

class FeffInpTemplate:
    def __init__(self, title, target_idx, cifloc) -> None:
        self.title = title
        self.target_idx = target_idx,
        self.cifloc = cifloc
        print(target_idx)
        self.string = [
            f"TITLE {self.title}",
            "",
            "EDGE K",
            "S02 1.0",
            "",
            "COREHOLE RPA",
            "XANES 4.0 0.07 0.0",
            "",
            "CONTROL 1 1 1 1 1 1",
            "PRINT 5 1 1 1 1 1",
            "",
            "SCF 4.0",
            "FMS 6.0",
            "*LDOS -30 15 0.01",
            "",
            "* Options for a k-space calculation :",
            "RECIPROCAL",
            "*REAL",
            "",
            "* Use 200 k-points:",
            "KMESH 200",
            "",
            f"TARGET {int(self.target_idx[0])}",
            "",
            "* This advanced option is not usually necessary:",
            "STRFAC 1.0 0.0 0.0",
            "",
            "* Read crystal structure from cif file:",
            f"CIF {self.cifloc}",
            "",
            "",
            "END",
            "",
        ]
        lines = "\n".join(self.string)
        with open("feff.inp", "w") as feff_inp:
            feff_inp.writelines(lines)