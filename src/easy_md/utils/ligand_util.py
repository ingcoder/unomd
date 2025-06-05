# Prepare ligand for forcefield parameterization


from openff.toolkit import Molecule
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff import ForceField


SUPPORTED = {"H","C","N","O","F","P","S","Cl","Br","I"}

def prepare_ligand_from_sdf(ligand_path: str) -> Molecule:
    lig = Molecule.from_file(ligand_path, allow_undefined_stereo=True)
    # lig = lig.with_aromaticity("MDL")

    # # element coverage check
    # if any(a.element.symbol not in SUPPORTED for a in lig.atoms):
    #     raise ValueError("Ligand contains elements outside Sage coverage")

    # ensure 3-D coords
    if not lig.conformers:
        try:
            lig.generate_conformers(n_conformers=1)
        except Exception as e:
            raise RuntimeError("RDKit failed to embed a conformer") from e

    # AM1-BCC charges
    try:
        lig.assign_partial_charges("am1bcc")
    except Exception as e:
        raise RuntimeError("AM1-BCC charge assignment failed") from e

    return lig
