import pandas as pd

df = pd.read_csv('file_name.csv')


def encode_atoms(triangles, method="old"):
    """
    Encode atoms as arrays representing the population in GA application.
    Write encoded items to files.

    :param triangles: List of triangles from CSV
    :param method: Name of method in which atoms are used, defines their encoding type
    """

    # TODO: urednost koda, OO?

    if method == "old":
        for triangle in triangles:
            genes = []

            # Nuc atom gene
            if triangle["Nuc_name"] == "OG":
                genes.append(0)
            elif triangle["Nuc_name"] == "SG":
                genes.append(1)

            # Acid atom gene
            if triangle["Acid_name"] == "OD1":
                genes.append(0)
            elif triangle["Acid_name"] == "OD2":
                genes.append(1)

            # Base atom gene
            if triangle["Base_aa"] == "HIS":
                genes.append(0)
            elif triangle["Base_aa"] == "ASP":
                genes.append(1)
            elif triangle["Base_aa"] == "GLU":
                genes.append(2)

            # - - - - - - - - - -

            # Nuc-Acid distance
            genes.append((1 - (8.57 - triangle["Dist_Nuc_Acid"]) / 2.84) * 10)

            # Acid-Base distance
            genes.append((1 - (4.62 - triangle["Dist_Acid_Base"]) / 1.96) * 10)

            # - - - - - - - - - -

            # Nuc angle
            genes.append((1 - (42.22 - triangle["Angle_Nuc"]) / 30.76) * 10)

            # Acid angle
            genes.append((1 - (65.85 - triangle["Angle_Acid"]) / 55.68) * 10)

            individual = ''.join(genes)
