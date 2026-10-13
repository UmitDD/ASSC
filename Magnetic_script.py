import ast

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
import spglib
import math

#--------------------------------IBZ folded data section--------------------
# IBZ folded data
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


def energy_fourier(basis,lat_cont,system,positions,types,mesh,scaling_factor):
    """

    :param basis: Number of magnetic basis interaction in the unit-cell. -Int-
    :param lat_cont: Lattice constant. -Float-
    :param system:Type of Bravais lattice (e.g: fcc, bcc). Type of Reciprocal lattice (e.g: fcc, bcc). -String-
    :param positions: Atoms' positions in the unit-cell
    :param types: Atoms' type in the unit-cell
    :param mesh: The unfolded mesh of propagation points (e.g: [10,10,10]) -List-
    :param scaling_factor: The scaling factor of Exchange energies (e.g: meV or mRy) -String-
    :return:
    """
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

        # =======Inter-atomic Distance======
        discase1 = pd.read_csv("dataoneone_alldirect.txt").sort_values(by=["hypo"]).reset_index()


        # ==== Set scaling factor based on desired output units ====#
        scaling = {"meV":1000,"mRyd":73.4985857}

        # ==== Tilted angle of magnetic moment along z-axis and magnetic moment of magnetic interaction basis ==== #

        mag_data = np.loadtxt("basis"+str(bs)+"/q0.00.00.0/mag_moments.csv",delimiter=",")
        mag_norm = np.linalg.norm(mag_data)
        mag_z = mag_data[-1]
        theta = np.degrees(np.arccos(mag_z/mag_norm))
        # print(f"Norm {mag_norm}")
        # print(f"mag_z {mag_z}")
        # print(f"theta {theta}")
        # ==== SS energies as function of propagation vectors ==== #
        energy_ss = data_unfold_ibz["energy"]
        energy_ss = energy_ss.values
        #print(f"ENERGY SS: {energy_ss}")

        # ==== Propagations vectors ==== #
        qx = data_unfold_ibz["qx"]; qy = data_unfold_ibz["qy"]; qz = data_unfold_ibz["qz"]
        q = np.array([qx, qy, qz]).T
        # ==== Inter-atomic distance in fractional coordinate ==== #
        R_atom = np.array(
            [np.round(discase1["dx"] / lat_cont, 2), np.round(discase1["dy"] / lat_cont, 2), np.round(discase1["dz"] / lat_cont, 2)])
        R_atom = R_atom.T
        # ==== Number of q-points and real-space points ==== #
        N_q = len(data_unfold_ibz["qx"])
        N_R = len(R_atom)
        # ==== Perform Fourier transformation to get Jij ==== #
        J_ij = np.zeros(N_R, dtype=np.complex128)
        J_R = []
        for exc_i, R in enumerate(R_atom):
            for qp, J_q_value in zip(q, energy_ss):
                # print(f"scaling {scaling[scaling_factor]}")
                # print(f"theta {np.sin(math.radians(theta))}")
                # print(f"energy[0] {energy_ss[0]}")
                # print(f"J_q {J_q_value}")
                phase_factor = np.exp(-1j * 2 * np.pi * np.dot(qp, R))
                J_ij[exc_i] += ((-1 * scaling[scaling_factor]) *  # H=-1*JiJj
                            (((J_q_value - energy_ss[0]) /
                              (np.sin(math.radians(theta)) * np.sin(math.radians(theta))))
                             * phase_factor))

            # ==== Normalize by the number of q-points ====
            J_ij[exc_i] /= N_q
            J_R.append([discase1["hypo"][exc_i], np.real(J_ij[exc_i])])
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
        dis = list(map(lambda num: num / lat_cont, dis))
        # ==== Calculate average J energies as terms of the interaction shell distance ==== #
        R = []
        for R_i in range(len(discase1["hypo"])):
            R.append([discase1["hypo"][R_i], 1])
        R_dict = {}
        for sublist in R:
            key = sublist[0]
            value = sublist[1]
            if key in R_dict:
                R_dict[key] += value
            else:
                R_dict[key] = value

        # ==== Print Exchange Values === #
        # for i, R in enumerate(R_atom):
        #     print(f"J(R = {R}) = {J_ij_real[i]:.6f}")

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
        # ax.set_yticks([-4, 0, 4, 8, 12, 16, 20, 24])  # Define specific y-axis tick positions
        # ax.set_xticks([0.5, 1, 1.5, 2, 2.5, 3, 3.5])  # Define specific y-axis tick positions
        # # Set limit
        # plt.xlim(0.4, 3.6)
        # plt.ylim(-5.005, 25)
        # # plt.ylim(-5.5, 25)
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
            s='bcc Fe',
            transform=ax.transAxes,
            fontsize=12,
            fontfamily='sans-serif',
            ha='left',
            va='top'
        )

        for line in ax.lines:
            line.set_clip_on(False)  # Disable clipping for markers at plot edges
        plt.savefig("output_basis"+str(bs)+"/"+"Exchange_graph.pdf", dpi=600, bbox_inches='tight', transparent=False)


        # ==== Curie Temperature Calculation =====#
        K_B = 0.0861733  # Kelvin over meV

        # ==== MFA- First equation --->KB*T_c = 2/3 * M**2 * Sum_j J_0j ==== #
        eq1 = (2 / 3) * sum(J_ij_real)
        Tc_MFA = eq1/K_B
        print(f"Tc_MFA {Tc_MFA}")
        # ==== RPA- Second equation ----> (T_c)^{-1} = 3KB/2N * sum sin^{2}(theta)/E(q)-E(0) ==== #
        eq2 = []
        for i in range(1, len(qx)):
            eq2.append((3 * K_B / (2 * len(qx))) *
                       ((np.sin(math.radians(theta) ** 2) /
                        (energy_ss[i] - energy_ss[0]))))

        Tc_RPA= 1 / sum(eq2)
        print(f"Tc_RPA {Tc_RPA}")
        data_Tc = pd.DataFrame([[Tc_MFA,Tc_RPA]], columns=["MFA", "RPA"])
        data_Tc.to_csv("output_basis" + str(bs) +"/"+"Tc_results", index=False)


if __name__ =="__main__":
    basis = int(sys.argv[1])
    lat_cont = float(sys.argv[2])
    system = str(sys.argv[3])
    positions = ast.literal_eval((sys.argv[4]))
    types = ast.literal_eval((sys.argv[5]))
    mesh = ast.literal_eval((sys.argv[6]))
    scaling_factor = str(sys.argv[7])

    energy_fourier(basis,lat_cont,system,positions, types, mesh, scaling_factor)

