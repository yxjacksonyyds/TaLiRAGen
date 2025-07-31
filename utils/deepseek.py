from openai import OpenAI

def chat(list):
    client = OpenAI(api_key='sk-0e18797e027e42e0bed0ff10c317940b', base_url='https://api.deepseek.com')
    dialogue_history = [
        {"role": "system", "content": "你是个乐于助人的分子生成专家。"}
    ]

    for i in list:
        user_input = i
        dialogue_history.append({"role": "user", "content": user_input})

        response = client.chat.completions.create(
            model="deepseek-reasoner",
            messages=dialogue_history,
            stream=False
        )
        assistant_message = response.choices[0].message.content
        dialogue_history.append({"role": "assistant", "content": assistant_message})
        # print("=================================================================")
        # print(f"Assistant: {assistant_message}")

    return assistant_message

def compare(answer_list,pref_name,binding_site):
    # Please install OpenAI SDK first: `pip3 install openai`
    client = OpenAI(api_key='sk-0e18797e027e42e0bed0ff10c317940b', base_url='https://api.deepseek.com')

    prompt = f'''作为拥有10年以上经验的药物设计专家，你需从针对{pref_name}靶点的候选分子中选出最优解。请系统执行以下分析流程：

**一、输入数据确认**
候选分子列表：{answer_list}
靶点特性：{pref_name}的蛋白结构特征、已知结合口袋性质、关键氨基酸残基信息{binding_site}

**二、多维评估体系**（权重分配）
1. 靶向特异性（40%）
   - 分子对接得分（使用AutoDock Vina标准）
   - 关键氢键形成能力
   - π-π堆积作用可能性
   - 疏水相互作用匹配度

2. ADMET特性（30%）
   - 类药五规则达标数
   - 血脑屏障渗透性预测
   - CYP450抑制概率
   - hERG毒性警示结构扫描

3. 合成路径可行性（20%）
   - 手性中心数量
   - 保护基需求度
   - 关键反应步骤收率预测
   - 原料商业可得性

4. IP潜力（10%）
   - SureChEMBL数据库新颖性检索
   - Markush结构覆盖分析

**三、决策树分析**
1. 初筛：排除具有PAINS结构或Rule of 3违规分子
2. 量化评分：对保留分子进行各维度0-10评分
3. 加权计算：总分=Σ(维度得分×权重)
4. 风险修正：对top3分子进行蒙特卡洛风险模拟

**四、输出要求**
▶ 最优分子：SMILES表达式（标注来源编号）
★ 优势矩阵：靶点结合能(ΔG)、QED评分、合成复杂度指数
⚠️ 风险预警：代谢不稳定位点标注、专利冲突提示
🔧 优化方案：建议引入的官能团及其位点
注意，一定要确保生成的分子smiles是合法的！！！！！
    '''
    
    response = client.chat.completions.create(
        model="deepseek-reasoner",
        messages=[
            {"role": "system", "content": prompt},
            {"role": "user", "content": "输出示例：<think>………………</think> <smiles>CCCC……</smiles>"},
        ],
        stream=False
    )
    
    print(response.choices[0].message.content)
    return response.choices[0].message.content

def chat_box(list):
    client = OpenAI(
      base_url="https://openrouter.ai/api/v1",
      api_key="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",
    )
    dialogue_history = [
        {"role": "system", "content": "你是个乐于助人的分子生成专家。"}
    ]

    for i in list:
        user_input = i
        dialogue_history.append({"role": "user", "content": user_input})

        response = client.chat.completions.create(
            model="deepseek/deepseek-r1:free",
            messages=dialogue_history,
            stream=False
        )

        assistant_message = response.choices[0].message.content
        dialogue_history.append({"role": "assistant", "content": assistant_message})
        # print("=================================================================")
        # print(f"Assistant: {assistant_message}")

    return assistant_message
