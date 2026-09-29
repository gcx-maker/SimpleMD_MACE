from mace.calculators import MACECalculator
from ase import Atoms


class MACEWrapper:

    def __init__(self, model_path):

        print("Loading MACE model...")

        self.calc = MACECalculator(
            model_paths=model_path,
            device="cpu",          # 如果用GPU改这里
            default_dtype="float32"
        )

        self.atoms = None

        print("MACE model loaded")


    def calculate(
        self,
        positions,
        cell,
        atomic_numbers
    ):

        # 第一次调用，创建ASE对象
        if self.atoms is None:

            self.atoms = Atoms(
                numbers=atomic_numbers,
                positions=positions,
                cell=cell,
                pbc=True
            )

            self.atoms.calc = self.calc


        else:

            # 后续只更新坐标
            self.atoms.positions[:] = positions

            # 如果你的MD是NPT，晶胞会变化
            # NVT可以不用这一行
            #self.atoms.cell[:] = cell


        energy = self.atoms.get_potential_energy()

        forces = self.atoms.get_forces()


        return energy, forces



# ===========================
# 全局只创建一次
# ===========================

mace = MACEWrapper(
    "IL_liquid_MACE_compiled.model"
)



def calculate(
    positions,
    cell,
    atomic_numbers
):

    return mace.calculate(
        positions,
        cell,
        atomic_numbers
    )