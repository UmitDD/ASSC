import ast
from vesin import NeighborList
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
import spglib
import math


def Lattice(lat_cont,system):
    """
    :param lat_cont: Lattice constant. -Float-
    :param Type of Bravais lattice (e.g: fcc, bcc).Type of Reciprocal lattice. Scaled to 2pi/lattice_constant(e.g: fcc_recip, bcc_recip) -String
    :return: Generate primitive unit-cell
    """
    dict = {"bcc":[[(lat_cont)*-0.5,(lat_cont)*0.5,(lat_cont)*0.5],
                [(lat_cont)*0.5,(lat_cont)*-0.5,(lat_cont)*0.5],
               [(lat_cont)*0.5,(lat_cont)*0.5,(lat_cont)*-0.5]],
            "fcc":[[0,(lat_cont)*0.5,(lat_cont)*0.5],
                    [(lat_cont)*0.5,0.0,(lat_cont)*0.5],
                    [(lat_cont)*0.5,(lat_cont)*0.5,0.0]],
            "bcc_recip":[[0.,  1.,  1.],[ 1., 0.,  1.], [ 1.,  1., 0.]],
            "fcc_recip":[[-1,  1.,  1.],[ 1., -1,  1.], [ 1.,  1., -1]]}


    return np.round(lat_cont, decimals=8), dict[system]

def interatomic(positions,unit_cell):

    calculator = NeighborList(cutoff=3,full_list=True)

    D,d = calculator.compute(points=positions,box=unit_cell,periodic=True,
                             quantities="Dd")


    distance_matrix = np.c_[(D, d)]
    return distance_matrix

def interatomic_two(lat_cont,system,delta):

    tol = 1e-8
    lattice_points = {"bcc":[np.array((0,0,0)),lat_cont*np.array((0.5,0.5,0.5))],
            "fcc":[np.array((0,0,0)),lat_cont*np.array((0.0,0.5,0.5)),
    lat_cont * np.array((0.5, 0.0, 0.5)),
    lat_cont * np.array((0.5, 0.5, 0.0))]}

    positions = lattice_points[system]
    unit_cell = lat_cont * np.eye(3)

    calculator = NeighborList(cutoff=lat_cont*3,full_list=True)

    i,j,D,d = calculator.compute(points=positions,box=unit_cell,periodic=True,
                             quantities="ijDd")


    mask = (i!=j) & (d>tol)
    R_frac = D[mask]/lat_cont
    R_frac = np.round(R_frac * 2.0)/2.0
    R_frac = np.unique(R_frac, axis=0)
    R_norm = np.linalg.norm(R_frac,axis=1)
    R_frac = R_frac[R_norm > tol]

    delta = np.asarray(delta,dtype=float)

    d_true = R_frac - delta
    true_dist = np.linalg.norm(d_true,axis=1)
    df = pd.DataFrame({
        "Rx": R_frac[:, 0],
        "Ry": R_frac[:, 1],
        "Rz": R_frac[:, 2],
        "true_dx": d_true[:, 0],
        "true_dy": d_true[:, 1],
        "true_dz": d_true[:, 2],
        "hypo": true_dist
    })

    df["R_norm"] = np.round(np.linalg.norm(R_frac,axis=1),6)
    df["True_norm"] = np.round(true_dist,6)
    df = df.sort_values(by=["R_norm","True_norm"]).reset_index(drop=True)

    nearest_true_by_R = (
        df
        .groupby("R_norm")["True_norm"]
        .min()
        .rename("nearest_true")
    )

    df = df.merge(nearest_true_by_R,on="R_norm",how="left")
    df_selected = df[df["True_norm"]==df["nearest_true"]].copy()
    df_selected = (
        df_selected
        .sort_values(by=["R_norm", "True_norm"])
        .reset_index(drop=True)
    )
    distance_matrix = np.c_[(df_selected[["Rx","Ry","Rz"]].values, df_selected["hypo"].values)]

    return distance_matrix


def read_poscar(filename):
    with open(filename) as f:
        lines = [ln.rstrip("\n") for ln in f]

    lattice_const = float(lines[1])

    unit_cell = np.array(
        [[float(x) for x in lines[i].split()] for i in range(2, 5)],
        dtype=float,
    )

    species = lines[5].split()
    counts  = [int(x) for x in lines[6].split()] # [2, 1, 1]
    n_atoms = sum(counts)

    positions = np.array(
        [[float(v) for v in lines[8 + i].split()] for i in range(n_atoms)],
        dtype=float,
    )

    return {
        "a": lattice_const,
        "unit_cell": unit_cell,
        "species": species,
        "counts": counts,
        "positions": positions,
    }

def exchange_fourier(basis,system,delta,types,mesh,scaling_factor,mag_basis,enSS90):
    """
    :param basis: Number of magnetic basis interaction in the unit-cell. -Int-
    :param system:Type of Bravais lattice (e.g: fcc, bcc). Type of Reciprocal lattice (e.g: fcc, bcc). -String-
    :param delta: Positional difference between two magnetic basis (atom) in the unit-cell.
    :param types: Atoms' type in the unit-cell
    :param mesh: The unfolded mesh of propagation points (e.g: [10,10,10]) -List-
    :param scaling_factor: The scaling factor of Exchange energies (e.g: meV or mRy) -String-
    :param mag_basis: Identify magnetic basis index in the POSCAR file
    :param enSS90: Spin-Spiral energy as function of phase 90 degrees between magnetic basis direction
    :return:
    """
    data_poscar = read_poscar("POSCAR")
    lat_cont = data_poscar["a"]
    positions = data_poscar["positions"]
    unit_cell = data_poscar["unit_cell"]


    for bs in range(1,basis+1):
        # ==== Read Energies of each magnetic interaction basis ==== #
        data = pd.read_csv("output_basis"+str(bs)+"/all_energies.txt", header=None, names=["qx", "qy", "qz","energy"])
        mask = data["energy"].isnull()
        if mask.any() == True:
            print("Energy data is missing. Check this directory(ies):\n")
            print(data[mask])
            sys.exit(1)
        else:
            print("All energies are present. The script is continuing")

        # ==== Create the cell and analyze symmetry ==== #
        lattice = Lattice(lat_cont, system)[1]
        cell = (lattice, positions, types)

        sym = spglib.get_symmetry_dataset(cell, symprec=1e-5)
        mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])

        # ==== Add Id column for each IBZ point ==== #
        data["id"] = np.unique(mapping)
        # ==== Create k-points and mapping to ir-grid points ==== #
        rec_lat = Lattice(lat_cont, system+"_recip")[1]
        wh_kp = []
        for index_i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
            wh_kp.append([ir_gp_id, np.round(np.dot(rec_lat, ((gp.astype(float) / mesh)[0]
                                                              , (gp.astype(float) / mesh)[1],
                                                              (gp.astype(float) / mesh)[2])), 2)])

        # ==== Unfold propagation points with their corresponding energies from IBZ representation ==== #
        unfold_ibz = []
        for i in range(len(wh_kp)):
            for j in range(len(data)):
                if wh_kp[i][0] == data["id"][j]:
                    unfold_ibz.append([wh_kp[i][1][0], wh_kp[i][1][1]
                                          , wh_kp[i][1][2], data["energy"][j]])
        # ==== Create dataset from unfolded propagation points and corresponding SS energies ==== #
        data_unfold_ibz = pd.DataFrame(unfold_ibz, columns=["qx", "qy", "qz", "energy"])
        data_unfold_ibz.to_csv("output_basis"+str(bs)+"/"+"energy_unfold", index=False)


        # ==== Set scaling factor based on desired output units ====#
        scaling = {"meV":1000,"mRyd":73.4985857}

        # ==== SS energies as function of propagation vectors ==== #
        energy_ss = data_unfold_ibz["energy"]
        energy_ss = energy_ss.values

        # ==== Propagations vectors ==== #
        qx = data_unfold_ibz["qx"];
        qy = data_unfold_ibz["qy"];
        qz = data_unfold_ibz["qz"]
        q = np.array([qx, qy, qz]).T
        if bs ==1 or bs ==2:
            # =======Inter-atomic Distance======
            discase1 = pd.DataFrame(interatomic(np.array(positions[mag_basis[0][bs - 1]]).reshape(1,-1), unit_cell),
                                    columns=["dx", "dy", "dz", "hypo"]).sort_values(by=["hypo"]).reset_index(drop=True)
            discase1.to_csv("interatomic_distance"+str(bs)+".csv", index=False)
            # ==== Inter-atomic distance in fractional coordinate ==== #
            R_atom = np.array(
                [np.round(discase1["dx"] / 1, 2), np.round(discase1["dy"] / 1, 2),
                 np.round(discase1["dz"] / 1, 2)])
            R_atom = R_atom.T
            # ==== Number of q-points and real-space points ==== #
            N_q = len(data_unfold_ibz["qx"])
            N_R = len(R_atom)

            # ==== Tilted angle of magnetic moment along z-axis and magnetic moment of magnetic interaction basis ==== #
            mag_data = pd.read_csv("basis" + str(bs) + "/q0.00.00.0/mag_moments.csv", sep=",", header=None)
            mag_data = mag_data.iloc[mag_basis[0][bs - 1]].values
            mag_norm = np.linalg.norm(mag_data)
            mag_z = mag_data[-1]
            theta = np.degrees(np.arccos(mag_z / mag_norm))
            np.savetxt("basis" + str(bs) + "/q0.00.00.0/theta.csv", [theta], delimiter=",")

            # ==== Perform Fourier transformation to get Jij ==== #
            J_ij = np.zeros(N_R, dtype=np.complex128)
            J_R = []
            for exc_i, R in enumerate(R_atom):
                for qp, J_q_value in zip(q, energy_ss):
                    phase_factor = np.exp(-1j * 2 * np.pi * np.dot(qp, R))
                    J_ij[exc_i] += ((-1 * scaling[scaling_factor]) *  # H=-1*JiJj
                                (((J_q_value - energy_ss[0]) /
                                  (np.sin(math.radians(theta)) * np.sin(math.radians(theta))))
                                 * phase_factor))

                # ==== Normalize by the number of q-points ====
                J_ij[exc_i] /= N_q
                J_R.append([np.round(discase1["hypo"][exc_i],4), np.round(np.real(J_ij[exc_i]),4)])

            # ==== Output real part of the exchange parameters J_ij ====
            J_ij_real = np.real(J_ij)

            J_dict = {}
            for sublist in J_R:
                key = sublist[0]
                value = sublist[1]
                if key in J_dict:
                    J_dict[key] += value
                else:
                    J_dict[key] = value

            # ==== Sort dictionary from the lowest inter-atomic distance to higher ==== #
            J_dict = dict(sorted(J_dict.items()))

            # Exchange energy as function of scaled distance
            dis = list(J_dict.keys())  # Shell radius
            # convert distance string list to float
            dis = list(map(float, dis))
            # Scale inter-atomic distance over lattice constant
            dis = list(map(lambda num: num / 1, dis))
            # ==== Calculate average J energies as terms of the interaction shell distance ==== #
            R = []
            for R_i in range(len(discase1["hypo"])):
                R.append([np.round(discase1["hypo"][R_i], 4), 1])
            R_dict = {}
            for sublist in R:
                key = sublist[0]
                value = sublist[1]
                if key in R_dict:
                    R_dict[key] += value
                else:
                    R_dict[key] = value

            # ==== Energy list in terms of Shell radius ====
            energy = []
            for key in R_dict:
                energy.append(J_dict[key] / R_dict[key])

            datatosave = pd.DataFrame({"distance": dis, "energy": energy})
            datatosave.to_csv("output_basis" + str(bs) + "/" + "Exchange_energies.csv", index=False)

            fig, ax = plt.subplots(figsize=(4.2, 4.2 / 1.618), dpi=600)  # Golden ratio
            ax.set_ylabel("$\mathrm{J_{ij}}$ (meV)", fontsize=12, fontname='sans-serif')
            ax.set_xlabel("R/a ($\mathrm{\AA}$)", fontsize=12, fontname='sans-serif')
            ax.tick_params(axis='both', which='major', length=2, direction='in')
            ax.axhline(y=0.0, color='lightgray', ls='-', linewidth=1)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
            ax.plot(dis, energy, "-", color="red", linewidth=1.5, marker='o', markersize=4,
                    mec='red', mfc='red', label="SS-VASP", zorder=11)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.legend(prop={'family': 'sans-serif', 'size': 9.5}, loc='upper right', frameon=True,
                       framealpha=1, edgecolor="none")

            # Add "bcc" label to the top-left corner (inside the plot)
            ax.text(
                x=0.2,
                y=0.95,
                s='Name',
                transform=ax.transAxes,
                fontsize=12,
                fontfamily='sans-serif',
                ha='left',
                va='top'
            )

            for line in ax.lines:
                line.set_clip_on(False)   # Disable clipping for markers at plot edges
            plt.savefig("output_basis" + str(bs) + "/" + "Exchange_graph.png", dpi=600, bbox_inches='tight',
                        transparent=False)
        if bs == 3:
            # =======Inter-atomic Distance======
            discase1 = pd.DataFrame(interatomic_two(lat_cont,system,delta),
                                    columns=["dx", "dy", "dz", "hypo"]).sort_values(by=["hypo"]).reset_index(drop=True)

            discase_new = discase1.drop_duplicates(["hypo"]).reset_index(drop=True)
            discase_new.to_csv("interatomic_distance" + str(bs) + ".csv", index=False)

            # ==== Inter-atomic distance in fractional coordinate ==== #
            R_atom = np.array(
                [np.round(discase1["dx"] / 1, 2), np.round(discase1["dy"] / 1, 2),
                 np.round(discase1["dz"] / 1, 2)])
            R_atom = R_atom.T

            # ==== Number of q-points and real-space points ==== #
            N_q = len(data_unfold_ibz["qx"])
            N_R = len(R_atom)

            enp90 = enSS90
            # ==== magnetic basis 1 ==== #
            basis_1 = pd.read_csv("output_basis1/energy_unfold")
            energy_ss1 =basis_1["energy"]
            mag_data = pd.read_csv("basis" + str(bs) + "/q0.00.00.0/mag_moments.csv", sep=",", header=None)

            mag_data1 = mag_data.iloc[mag_basis[0][0]].values
            mag_norm1 = np.linalg.norm(mag_data1)
            mag_z1 = mag_data1[-1]
            theta1 = np.degrees(np.arccos(mag_z1 / mag_norm1))

            basis_2 = pd.read_csv("output_basis2/energy_unfold")
            energy_ss2 = basis_2["energy"]

            mag_data2 = mag_data.iloc[mag_basis[0][1]].values
            mag_norm2 = np.linalg.norm(mag_data2)
            mag_z2 = mag_data2[-1]
            theta2 = np.degrees(np.arccos(mag_z2 / mag_norm2))

            # ==== Perform Fourier transformation to get Jij ==== #
            J_ij = np.zeros(N_R, dtype=np.complex128)
            J_R = []
            for exc_i, R in enumerate(R_atom):
                for qp, J_q12_value, J_q11_value, J_q22_value in zip(q, energy_ss,
                                                                     energy_ss1, energy_ss2):
                    phase_factor = np.exp(-1j * 2 * np.pi * np.dot(qp, R))
                    J_ij[exc_i] += ((0.5* scaling[scaling_factor]) *  # H=-1*JiJj
                                ((((enp90 - J_q12_value) / (np.sin(math.radians(theta1)) * np.sin(math.radians(theta2)))
                                   + 1 * (J_q11_value - energy_ss1[0]) / (np.sin(math.radians(theta1)) * np.sin(math.radians(theta2)))
                                   + 1 * (J_q22_value - energy_ss2[0]) / (np.sin(math.radians(theta1)) * np.sin(math.radians(theta2)))))
                                 * phase_factor))
                # ==== Normalize by the number of q-points ====
                J_ij[exc_i] /= N_q
                J_R.append([np.round(discase1["hypo"][exc_i], 4), np.round(np.real(J_ij[exc_i]), 4)])

            # ==== Output real part of the exchange parameters J_ij ====
            J_ij_real = np.real(J_ij)

            J_dict = {}
            for sublist in J_R:
                key = sublist[0]
                value = sublist[1]
                if key in J_dict:
                    J_dict[key] += value
                else:
                    J_dict[key] = value

            # ==== Sort dictionary from the lowest inter-atomic distance to higher ==== #
            J_dict = dict(sorted(J_dict.items()))

            # Exchange energy as function of scaled distance
            dis = list(J_dict.keys())  # Shell radius
            # convert distance string list to float
            dis = list(map(float, dis))
            # Scale inter-atomic distance over lattice constant
            dis = list(map(lambda num: num / 1, dis))
            # ==== Calculate average J energies as terms of the interaction shell distance ==== #
            R = []
            for R_i in range(len(discase1["hypo"])):
                R.append([np.round(discase1["hypo"][R_i],4), 1])
            R_dict = {}
            for sublist in R:
                key = sublist[0]
                value = sublist[1]
                if key in R_dict:
                    R_dict[key] += value
                else:
                    R_dict[key] = value

            # ==== Energy list in terms of Shell radius ====
            energy = []
            for key in R_dict:
                energy.append(J_dict[key] / R_dict[key])

            datatosave = pd.DataFrame({"distance": dis, "energy": energy})
            datatosave.to_csv("output_basis"+str(bs)+"/"+"Exchange_energies.csv", index=False)

            fig, ax = plt.subplots(figsize=(4.2, 4.2 / 1.618), dpi=600)  # Golden ratio
            ax.set_ylabel("$\mathrm{J_{ij}}$ (meV)", fontsize=12, fontname='sans-serif')
            ax.set_xlabel("R/a ($\mathrm{\AA}$)", fontsize=12, fontname='sans-serif')
            ax.tick_params(axis='both', which='major', length=2, direction='in')
            ax.axhline(y=0.0, color='lightgray', ls='-', linewidth=1)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
            ax.plot(dis, energy, "-", color="red", linewidth=1.5, marker='o', markersize=4,
                    mec='red', mfc='red', label="SS-VASP", zorder=11)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.legend(prop={'family': 'sans-serif', 'size': 9.5}, loc='upper right', frameon=True,
                       framealpha=1, edgecolor="none")

            # # Add "bcc" label to the top-left corner (inside the plot)
            # ax.text(
            #     x=0.2,
            #     y=0.95,
            #     s='Name',
            #     transform=ax.transAxes,
            #     fontsize=12,
            #     fontfamily='sans-serif',
            #     ha='left',
            #     va='top'
            # )

            for line in ax.lines:
                line.set_clip_on(False)  # Disable clipping for markers at plot edges
            plt.savefig("output_basis"+str(bs)+"/"+"Exchange_graph.png", dpi=600, bbox_inches='tight', transparent=False)


        # ==== Curie Temperature Calculation =====#
        # K_B = 0.0861733  # Kelvin over meV
        #
        # # ==== MFA- First equation --->KB*T_c = 2/3 * M**2 * Sum_j J_0j ==== #
        # eq1 = (2 / 3) * sum(J_ij_real)
        # Tc_MFA = eq1/K_B
        # print(f"Tc_MFA {Tc_MFA}")
        # # ==== RPA- Second equation ----> (T_c)^{-1} = 3KB/2N * sum sin^{2}(theta)/E(q)-E(0) ==== #
        # eq2 = []
        # for i in range(1, len(qx)):
        #     eq2.append((3 * K_B / (2 * len(qx))) *
        #                ((np.sin(math.radians(theta) ** 2) /
        #                 (energy_ss[i] - energy_ss[0]))))
        #
        # Tc_RPA= (1 / sum(eq2)) * 1000
        # print(f"Tc_RPA {Tc_RPA}")
        # data_Tc = pd.DataFrame([[Tc_MFA,Tc_RPA]], columns=["MFA", "RPA"])
        # data_Tc.to_csv("output_basis" + str(bs) +"/"+"Tc_results", index=False)


if __name__ =="__main__":
    basis = int(sys.argv[1])
    system = str(sys.argv[2])
    delta = ast.literal_eval(sys.argv[3])
    types = ast.literal_eval((sys.argv[4]))
    mesh = ast.literal_eval((sys.argv[5]))
    scaling_factor = str(sys.argv[6])
    mag_basis = ast.literal_eval(sys.argv[7])
    enSS90 = float(sys.argv[8])
    exchange_fourier(basis,system,delta,types, mesh, scaling_factor,mag_basis,enSS90)

