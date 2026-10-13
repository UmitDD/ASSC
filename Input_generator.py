import numpy as np
import pandas as pd
import ast
import sys
import os
import scipy.linalg._solvers as la

def bravais_lattice(lat_cont,system):
    """

    :param lat_cont: Lattice constant. -Float-
    :param system: Type of Bravais lattice (e.g: fcc, bcc). -String-
    :return: Generate primitive unit-cell
    """
    dict = {"bcc":[[(lat_cont)*-0.5,(lat_cont)*0.5,(lat_cont)*0.5],
                [(lat_cont)*0.5,(lat_cont)*-0.5,(lat_cont)*0.5],
               [(lat_cont)*0.5,(lat_cont)*0.5,(lat_cont)*-0.5]],
            "fcc":[[0,(lat_cont)*0.5,(lat_cont)*0.5],
                    [(lat_cont)*0.5,0.0,(lat_cont)*0.5],
                    [(lat_cont)*0.5,(lat_cont)*0.5,0.0]]}



    return np.round(lat_cont, decimals=8), dict[system]



def rec(ai, aj, ak):
    """
    Reciprocal lattice generator
    :param ai: lattice constant at x-coordinate. -Float-
    :param aj: lattice constant at y-coordinate. -Float-
    :param ak: lattice constant at z-coordinate. -Float-
    :return: Reciprocal lattice along x,y, and z direction in Brillouin zone
    """
    b = 2 * np.pi * (np.cross(aj, ak) / (np.dot(ai, np.cross(aj, ak))))
    return b

def incar_generate(basis,pro_vectors,lat_cont,system,mag_moment,const_moment,file_name):
    """

    :param basis: Number of magnetic basis interaction in the unit-cell. -Int-
    :param pro_vectors:Input file of Propagation path has three column: qx,qy,qz. The coordinate representation of
    propagation vectors are fractional reciprocal lattice of primitive unit-cell. -String-
    :param lat_cont: Lattice constant. -Float-
    :param system: Type of Bravais lattice (e.g: fcc, bcc). -String-
    :param mag_moment:
    :param const_moment:
    :param file_name:
    :return:
    """
    # ==== Propagation Vectors Input Files ==== #
    data_pro = pd.read_csv(pro_vectors, sep=",")
    qx = data_pro.iloc[:, 0].values.tolist()
    qy = data_pro.iloc[:, 1].values.tolist()
    qz = data_pro.iloc[:, 2].values.tolist()


    lattice_system = bravais_lattice(lat_cont, system)
    lattice = lattice_system[0]
    bravais = lattice_system[1]
    # ==== Reciprocal lattice matrix in fractional coordinate ==== #
    rp_m = [(lattice / (2 * np.pi)) * rec(bravais[0], bravais[1], bravais[2]), (lattice / (2 * np.pi)) * rec(bravais[1], bravais[2], bravais[0]),
            (lattice / (2 * np.pi)) * rec(bravais[2], bravais[0], bravais[1])]

    # ==== Solve linear system for fractional representation of propagation vector in terms of unit-cell symmetry ==== #
    qf = []
    for i in range(len(qx)):
        qf.append(la.solve([rp_m[0], rp_m[1], rp_m[2]], [qx[i], qy[i], qz[i]]))
    # ==== Create csv file from scaled q-parameters to primitive unit cell coordinate system in reciprocal space ==== #
    data = pd.DataFrame(np.round(qf, 6), columns=["qx", "qy", "qz"])
    data.to_csv("propagation_vectors_irreducible_primitive.csv", index=False)

    # ==== Write INCAR Files ==== #
    with open("example_INCAR","r") as file:
        content = file.read()

    # ==== Write Job Script to execute SS calculation ==== #
    with open("example_scr.sh","r") as scr:
        content_scr = scr.read()

    # ==== Write Cat Script to merge SS energies into format that use in exchange calculation==== #
    with open("example_cat.sh", "r") as cat:
        content_cat = cat.read()

    for at in range(1,basis+1):
        path = './INCAR_files_basis' + str(at)
        os.makedirs(path, exist_ok=True)

        const = " ".join(str(c) for con in const_moment[at-1] for c in con)
        mag = " ".join(str(m) for magm in mag_moment[at-1] for m in magm)

        lines = [f"M_CONSTR={const}\n",
                f"MAGMOM={mag}\n" ]

        # ==== QSPIRAL tags for each magnetic basis ==== #
        for i in range(0,len(qx)):

            qsprial_tag=("\nQSPIRAL =" + str(np.round(qf[i][0], 6)) + " "
                         + str(np.round(qf[i][1], 6)) + " "
                         + str(np.round(qf[i][2], 6)) + "\n")

            new_content = content + "\n"+"".join(lines) + qsprial_tag

            file_path = (str(path)+"/"+file_name+"q"+str(qx[i]) + str(qy[i]) + str(qz[i]))

            with open(file_path, "w") as f:
                f.write(new_content)

            # ==== Spin Spiral Job Calculation Script ==== #
            scr_lines = ["qx=(" + ' '.join(["{:}".format(value) for value in qx]) + ")" + "\n",
                 "qy=(" + ' '.join(["{:}".format(value) for value in qy]) + ")" + "\n",
                 "qz=(" + ' '.join(["{:}".format(value) for value in qz]) + ")" + "\n",
                "mkdir basis" + str(at) + "\n"
                , "cd basis" + str(at) + "\n"
                , "for index in" + " " + '"${!qy[@]}"' + "\n"
                , "do" + "\n"
                , 'k="${qx[index]}"' + "\n"
                , 'i="${qy[index]}"' + "\n"
                , 'j="${qz[index]}"' + "\n"
                , 'dir_name="q${k}${i}${j}"' + "\n"
                , 'if [ -d "$dir_name" ]; then' + "\n"
                , 'echo "$dir_name already exists... skipping"' + "\n"
                , 'else' + "\n"
                , 'mkdir "$dir_name"' + "\n"
                , 'cd "$dir_name" || exit' + "\n"
                , "   cp ../../INCAR_files_basis" + str(at) + "/"+file_name + '"$dir_name"' + ' INCAR\n'
                , "   cp ../../KPOINTS ." + "\n"
                , "   cp ../../POTCAR ." + "\n"
                , "   cp ../../POSCAR ." + "\n"
                , "   mpirun $VASP   >> vasp.out 2> err" + "\n"
                , "   en=$(grep 'TOTEN' OUTCAR | tail -n 1 | awk '{print $5}')\n"
                , '   echo "${k},${i},${j},$en" > SS_${k}${i}${j}.txt\n'
                , "   [ -f vasprun.xml ] && grep 'QSPIRAL' vasprun.xml >> qspiral || true\n"
                , "   last_hdr=$(grep -n 'ion[[:space:]]\\+MW_int[[:space:]]\\+M_int' OSZICAR | tail -1 | cut -d: -f1)\n"
                , "   [ -z \"$last_hdr\" ] && { echo 'No M_int block in OSZICAR' >&2; exit 1; }\n"
                , r'''awk -v s="$last_hdr" 'NR>s { if ($1 ~ /^[0-9]+$/ && $2 ~ /^-?[0-9.]+$/ && NF>=7) printf "%s,%s,%s\n", $5, $6, $7; else exit }' OSZICAR > mag_moments.csv'''+"\n"
                , "   rm vasprun.xml" + "\n"
                , "   rm CHGCAR" + "\n"
                , "   rm CHG" + "\n"
                , "   rm PROCAR" + "\n"
                , "   rm POTCAR" + "\n"
                , "   rm DOSCAR" + "\n"
                , "   cd ../ || exit" + "\n"
                , "fi" + "\n"
                , "done" + "\n"]

            # ==== Energy merge (cat) Script ==== #

            cat_lines = [f"mkdir -p output_basis{at}\n",
                    "qx=(" + ' '.join(["{:}".format(value) for value in qx]) + ")" + "\n",
                     "qy=(" + ' '.join(["{:}".format(value) for value in qy]) + ")" + "\n",
                     "qz=(" + ' '.join(["{:}".format(value) for value in qz]) + ")" + "\n"
                ,"for index in" + " " + '"${!qy[@]}"' + "\n"
                , "do" + "\n"
                , ' k="${qx[index]}"' + "\n"
                , ' i="${qy[index]}"' + "\n"
                , ' j="${qz[index]}"' + "\n"
                , f' cp basis{at}/'+'q${k}${i}${j}'+'/SS_${k}${i}${j}.txt'+f' output_basis{at}'"\n"
                , f' file="output_basis{at}'+'/SS_${k}${i}${j}.txt"'+"\n"
                ,   'if [ -f "$file" ]; then' + "\n"
                ,   'cat "$file" >>'+f'output_basis{at}/'+'all_energies.txt'+ "\n"
                ,   'else' + "\n"
                ,   'echo SS_${k}${i}${j}.txt is missing >>'+f'output_basis{at}/'+'missing_energy_files.txt\n'
                ,   'fi\n'
                , 'done'
                ]

            new_content_scr = content_scr + "\n" + "".join(scr_lines)
            with open("job_SS"+str(at)+".sh", "w") as f:
                f.write(new_content_scr)

            new_content_cat = content_cat + "\n" + "".join(cat_lines)
            with open("cat" + str(at)+".sh", "w") as f:
                f.write(new_content_cat)



if __name__=="__main__":
    basis = int(sys.argv[1])
    pro_vectors = sys.argv[2]
    lat_cont = float(sys.argv[3])
    system = str(sys.argv[4])
    mag_moment = ast.literal_eval(sys.argv[5])
    const_moment = ast.literal_eval(sys.argv[6])
    file_name = sys.argv[7]

    incar_generate(basis, pro_vectors, lat_cont, system, mag_moment, const_moment, file_name)