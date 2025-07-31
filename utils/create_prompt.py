class promptt:
    def __init__(self,pref_name,binding_site,Ki,request,chemblid,text):
        self.pref_name = pref_name
        self.binding_site = binding_site
        self.Ki = Ki
        self.request = request
        self.chemblid = chemblid
        self.text = text

    def rule(self):
        prompt = '''SMILES生成与校验结构化指南
（分步骤任务执行模板）

一、原子规范

元素符号规则

单原子：C O N → 正常价态无需括号

双原子：Cl Br → 第二字母必须小写

非常见价态：[Fe+3] [S-2]

芳香原子：c1ccccc1（苯环）

!错误检测!
× Cl → ✗
√ Cl → ✓

二、化学键处理流程

键类型优先级
默认单键：CC（乙烷）
双键：C=C（乙烯）
三键：C#C（乙炔）
芳香键：c:c（隐式）

立体化学标记
顺式：F/C=C/F
反式：F\C=C\F
四面体手性：C@Br → 顺时针
C@@Br → 逆时针

三、结构构建步骤

主链构建
线性结构：CCCC（丁烷）
环状结构：C1CCCC1（环戊烷）

分支处理
单分支：CC(C)C（异丁烷）
多分支：CC(C)(C)C（新戊烷）

四、特殊元素处理规范

氢原子标注
必须显式的情况：

带电：[H+]

同位素：[2H]

单质：[H][H]

电荷表示标准
单电荷：[Na+]
多电荷：[Fe+3] 或 [Fe+++]

五、复合结构处理

分离组分：用.分隔
[Na+].[Cl-]（氯化钠）
CC(=O)O.CCO（乙酸+乙醇）

环编号规则

成对数字标记：C1CCCCC1

避免数字冲突：C1CC2CCCC2C1

六、验证检查清单

语法校验点
□ 所有非标准价态元素是否加[]
□ 芳香环是否闭合（如c1ccccc1）
□ 环编号是否成对出现
□ 立体标记是否配套（@/@@，/与\）

语义校验点
□ 氢原子数量是否符合价键规则
□ 电荷是否平衡（整体中性优先）
□ 同位素标记是否正确（如[13C]）

七、示例生成模板
输入结构：环己烷（椅式构型）
正确SMILES：C1CCCCC1

输入结构：顺式二氯乙烯
正确SMILES：Cl/C=C/Cl

输入结构：硫酸根离子
正确SMILES：S(=O)([O-])[O-]
        '''
        return prompt
        
    def get_binding_site_propmt(self):
        prompt = f'''###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. Now, you have accepted the literature and summarized some knowledge from it. You can obtain the knowledge you have learned through context. Your current stage: binding site analysis (structural docking), studying the binding site information of the target I provided to you.
###Role: Molecular Generation Expert (Structural Analysis Mode)
###Specific task: At this stage, your task is to analyze some binding site related information of the target downloaded from the uniprotoKB database in JSON format, learn the structure related knowledge of the target from you, and use it as the context for subsequent molecular generation.
###Input: {self.binding_site} (UniProt exported binding site data)
###Analysis steps:
1. Extract structural features that combine positions
2. Identify the molecules or residues that can bind to the binding site and analyze their characteristics
3. Try to piece together the molecular residues bound to that position into a single molecule
4. Macroscopically provide the spatial significance of the binding site of the target, laying the foundation for the subsequent generation of bound molecules.
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, in order to better learn the knowledge of the binding site of the target mentioned above.)
###Output: Knowledge learned
'''
        return prompt


    def get_resorted_prompt(self):
        # 将文本内容放在json材料后
        prompt1 = self.get_article()
        prompt2 = self.get_binding_site_propmt()
        prompt3 = self.get_Ki_ligand()
        prompt4 = self.get_generate_prompt()
        prompt5 = self.get_check_prompt()
        rule = self.rule()
        return rule,prompt2,prompt3,prompt1,prompt4,prompt5,self.chemblid

    def get_icnm(self,icnm):
        prompt = f'''###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. You can obtain the knowledge you have learned through context. Your current stage: Experimental validation has been combined with ligand analysis to study the information on the existing binding ligands of the target that I have provided to you.
###Role: Molecular Generation Expert (Structural Analysis Mode)
###Specific task: At this stage, your task is to analyze some highly active ligand related information of the target downloaded from the CHEMBL database in JSON format, learn about the ligand related knowledge of the target from you, and use it as the context for subsequent molecule generation.
###Input: {icnm} (icnm binding site data exported from UniProt)
###Analysis steps:
1. Establish a 2D-QSAR model: identify key descriptors that affect pIC50 (such as AlogP, TPSA)
2. Conduct 3D pharmacophore comparison: Calculate molecular overlap RMSD and extract conservative features
3. Building an active cliff analysis: identifying substituent patterns that lead to sudden changes in activity
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, with the aim of better learning about the relevant ligands of the target mentioned above.)
###Output: Knowledge learned
'''
        return prompt

    def get_icnm_prompt(self,icnm):
        prompt1 = self.get_article()
        prompt2 = self.get_binding_site_propmt()
        prompt3 = self.get_icnm(icnm)
        prompt4 = self.get_generate_prompt()
        prompt5 = self.get_check_prompt()
        rule = self.rule()
        return rule,prompt1,prompt2,prompt3,prompt4,prompt5,self.chemblid

    def get_article(self):
        prompt = f'''###Background: You are an expert in the field of molecular generation and can proficiently generate ligand molecules based on some information of proteins. Smiles said that you will gain relevant knowledge step by step through multiple rounds of conversations with me, and gradually improve your ability to generate target ligand molecules. Your current stage: literature knowledge structuring (knowledge anchoring), you can conduct ligand related research on a certain protein target based on the following literature information about the target.
###Role: Molecular Generation Expert (Literature Analysis Mode)
###Input: {self.text} (target information from paper/network)
###Task: Perform the following analysis process:
1. Basic information analysis of the target
2. Analysis of structural information of the target
3. Analysis of experimental information on the known target
4. Analysis of ligand generation mode of the target
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, with the aim of better learning the knowledge of the target article mentioned above.)
###Output: Knowledge obtained from the article.
        '''
        return prompt
    
    def get_Ki_ligand(self):
        prompt = f'''
###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. You can generate the smiles representation of the ligand molecule for a certain protein target based on the following literature information, binding sites, existing drug molecules, and active substances with high binding affinity. Your current stage: ligand SAR analysis (structure-activity relationship modeling)
###Task: At this stage, your task is to analyze the ligand related information of some high binding affinity values (arranged in descending order of binding ability) of the target downloaded from the Chembl database in JSON format, and learn the relevant knowledge of the target from you as the context for subsequent molecule generation.
###Role: Molecular Generation Expert (QSAR Analysis Mode)
###Input: {self.Ki} (ChEMBL High Activity Molecular Dataset)
###Analysis process:
1. Extract common features of ligands
2. Analyze the binding mode between ligands and targets
3. Extract the molecular fragments that make up the ligand in the example
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, with the aim of better learning about the relevant ligands of the target mentioned above.)
###Output: Knowledge learned
        '''
        return prompt
    
    def get_generate_prompt(self):
        prompt = f'''###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. You can generate the smiles representation of the ligand molecule for a certain protein target based on the following literature information, binding sites, existing drug molecules, and active substances with high binding affinity. Current stage: Molecular generation (rational design), designing ligand molecules for the target.
###Task: At this stage, your task is to generate new ligand molecules for the target based on the knowledge you have obtained, in order to obtain molecules with higher binding affinity, higher drug like properties, and higher molecular synthesizability compared to the provided knowledge. Please generate the ligand molecule for the target based on the knowledge of "characteristics of high affinity molecules" and "methods for generating high affinity molecules smiles" in the context, think step by step, provide the reasoning process, check the reasoning process, and summarize the final smile results. Require an affinity Vina score below -9.
###Role: Molecular Generation Expert (de novo design pattern)
###Molecular design logic:
<Step 1>Find the core structure in the ligand as the reference molecule.
<Step 2>Introduce drug-resistant groups (such as rigid macrocyclic structures).
<Step 3>Optimize the molecule for easier synthesis.
Ensure the generation of legal and high binding molecules.
###Output specification:
xml
<design-process>
<base>Reference molecule</base>
<step1>
<modify>First modified content</modify>
<smiles>First modified smiles</smiles>
</step1>
<step2>
<modify>Second modification content</modify>
<smiles>Second modified smiles</smiles>
</step2>
<final>Final smiles</final>
</design-process>
###注意：请你生成最大环尺寸在8以下的分子。
        '''
        return prompt
    
    def get_check_prompt(self):
        prompt = f'''Background: You are an expert in the field of molecular generation and can proficiently generate ligand molecules for proteins based on certain information. You can generate the smiles representation of the ligand molecule for a certain protein target based on the following literature information, binding sites, existing drug molecules, and active substances with high binding affinity. Your current stage is molecular validation (multidimensional evaluation), where you need to evaluate the molecules generated in the previous section, validate, rank, and output them.
###Task: At this stage, your task is to carefully read the molecules you generate in the context, perform verification, validation, and ranking work, and ensure that the smiles are valid. It is forbidden to generate invalid molecules!!! Definitely generate molecules with strong binding affinity!!!! Among them, the molecule ranked first must be the one with the strongest binding force, and we must spare no effort to improve its binding ability, while maintaining its legitimacy.
###Verification Protocol:
1. Chemical validity check: Use RDKit to verify the correctness of SMILES syntax/valence state.
2. Docking verification: AutoDock Vina rapid docking with the original crystal structure (Δ G ≤ -9 kcal/mol)
3. Drug Evaluation: QED ≥ 0.6 and SAscore ≤ 4
4. Comprehensive sorting: Weighted scores based on (0.4 × pKiupred)+(0.3 × QED)+(0.3 × SAscore)
###Output template:
xml
<validation-summary>
<smiles>The only legitimate and highly binding molecule for testing</smiles>
<qed>Predict its qed indicators</qed>
<affinity>Predict their Vina score</affinity>
<sa>Predict its synthesizability index sa</sa>
</validation-summary>
###Attention: The entire process must ensure the generation of legal smiles!!!!!!!!!! The conditions for legal smiles:
1. Grammar level: Atomic symbols must be valid (such as C, [Na+]), bond connections must be clear (single bonds can be omitted, double/triple bonds use=#), parentheses are closed in pairs, ring markers match numbers (such as C1CCCC1), aromatic atoms are lowercase (such as benzene c1ccccc1); 2. At the chemical level: the atomic valence state is reasonable (such as carbon not exceeding 4 bonds), the hydrogen atom is implicitly expressed (special explicit hydrogen such as [NH3]), the charge labeling is correct (such as [Fe+2]), and the overall structure is free of contradictions (such as ring energy closure and ion charge balance).
###注意：请尽量避免以下情况的分子生成：①.出现分子数超过8的大环 ②.出现三个环共用一个分子的情况 ③.分子中杂原子太多 ④.分子中杂原子在同一个环中多次出现。
'''
        return prompt
    
    def get_prompt(self):
        prompt1 = self.get_article()
        prompt2 = self.get_binding_site_propmt()
        prompt3 = self.get_Ki_ligand()
        prompt4 = self.get_generate_prompt()
        prompt5 = self.get_check_prompt()
        # smiles = self.rule()
        return prompt1,prompt2,prompt3,prompt4,prompt5,self.chemblid

    def get_prompt_from_id(self,uniprotid,chemblid):
        # 设置 URL
        url = f"https://www.ebi.ac.uk/chembl/api/data/target/{chemblid}.json"
        
        # 发送 GET 请求
        response = requests.get(url)
        
        # 检查请求是否成功
        if response.status_code == 200:
            # 解析 JSON 数据
            data = response.json()
            return data['pref_name']
        else:
            return None

    def get_function(self,function,paper):
        prompt = f'''你好，你是一个经验丰富的蛋白质靶点的配体生成专家。现在我会交给你{self.pref_name}蛋白质靶点的功能信息，请你研究这些信息：{function}。在研究完结构信息后，请你阅读目前与之相关的文献：{paper}。好了，在阅读完文献并体会过该靶点的功能信息后，请你尝试生成该靶点结合的配体的分子smiles表示。最终结果使用<smiles></smiles>包裹，以便后续提取。
        '''

    def get_youhua(self,smiles,score,constrain):
        prompt = '''Role: You are an expert computational medicinal chemist specializing in ligand optimization. Your task is to analyze existing ligand candidates and generate improved novel molecules by intelligently combining their strengths.
Input Data Structure:
molecules: [list of SMILES strings]
overall_scores: [list of integers 0-100] ← Higher is better for drug-likeness
constraint_metrics: {
"metric_name": {
"values": [list of numbers],
"direction": "lower" or "higher" // e.g., "ring_count ≤5" → direction="lower"
},
...
}
Core Instructions:
Analyze molecular strengths
Identify top 3 ligands with the highest overall scores and extract their drug-like features (e.g., hydrogen bond donors, aromatic rings)
For each constraint metric:
Find ligands optimally satisfying direction (e.g., lowest values if direction="lower")
Extract structural elements enabling constraint compliance
Cross-reference to discover:
Why high-scoring ligands failed constraints
Why constraint-compliant ligands scored poorly overall
Generate novel ligands
Combine extracted features to create 3-5 new SMILES strings that:
Integrate complementary strengths from multiple input ligands
Explicitly optimize for both high overall scores AND constraint satisfaction
For each new ligand:
Explain which input ligand features were fused (e.g., "torsional flexibility from Ligand#3 + polar surface area from Ligand#1")
Predict improvements:
overall_score rationale (e.g., "Added sulfonyl group to enhance solubility")
constraint rationale (e.g., "Reduced rotatable bonds to meet direction='lower' requirement")
Output Requirements:
### Optimized Ligand: <SMILES>The refined ligand molecule</SMILES>  
        '''
    def get_no_rag(self):
        prompt = f'''请你生成{self.pref_name}靶点的配体分子smiles，使用<smiles></smiles>包裹，以便后续提取'''
        return prompt

        
