from utils.evaluation.docking_vina import PrepLig,PrepProt,VinaDock
from rdkit import Chem
from rdkit.Chem import QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Descriptors
import math

class test:
    def __init__(self):
        pass  
    def cauculate(self,smiles,chemblid):
        lig_pdbqt = '/root/IRDiff/utils/evaluation/data/ligand.pdbqt'
        mol_file = smiles
        a = PrepLig(mol_file, 'smi')
        a.addH()
        a.gen_conf()
        a.get_pdbqt(lig_pdbqt)
        
        prot_file = f'/root/autodl-tmp/pdb_files/{chemblid}.pdb'
        prot_dry = '/root/IRDiff/pre_deal_for_generate/data/protein_dry.pdb'
        prot_pqr = '/root/IRDiff/pre_deal_for_generate/data/protein.pqr'
        prot_pdbqt = '/root/IRDiff/pre_deal_for_generate/data/protein.pdbqt'
        b = PrepProt(prot_file)
        b.del_water(prot_dry)
        b.addH(prot_pqr)
        b.get_pdbqt(prot_pdbqt)
        
        dock = VinaDock(lig_pdbqt, prot_pdbqt)
        dock.get_box()
        result = dock.dock()
        print(f"vina分数是{result}")
        # 计算指标
        qed_value = self.calculate_qed(mol_file)
        sa_value = self.calculate_sa(mol_file)
        
        return qed_value,sa_value,result
        
        
       
        
    def calculate_qed(self,smiles):
        molecule = Chem.MolFromSmiles(smiles)
        return QED.qed(molecule)
    
    def calculate_sa(self,smiles):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError("Invalid SMILES string")
        
        # 使用分子特征来估算 SA
        num_atoms = molecule.GetNumAtoms()
        num_rings = rdMolDescriptors.CalcNumRings(molecule)
        num_chiral_centers = len(Chem.FindMolChiralCenters(molecule, includeUnassigned=True))
        
        # 简单的 SA 评分公式
        sa_score = 1.0 + math.log(num_atoms + num_rings + num_chiral_centers + 1)
        return sa_score
    
    



