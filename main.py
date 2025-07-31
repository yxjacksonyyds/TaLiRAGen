import json
from utils.get_knowledge import get_bind,get_activity
from utils.create_prompt import promptt
import requests
from utils.qwen import get_answer
from transformers import pipeline
from modelscope import AutoModelForCausalLM, AutoTokenizer
from accelerate import init_empty_weights  # 新路径

with open('dataset/filtered_TTD.txt', 'r') as data:
    uniprot_list = []
    chembl_list = []
    for i in data:
        list = i.split(' ')
        uniprot_list.append(list[0])
        chembl_list.append(list[1])
# 读取数据集，并加载ID

with open("data_process/process.txt", "w") as file:
    for uniprotid, chemblid in zip(uniprot_list, chembl_list):
        if get_bind(uniprotid) and get_activity(chemblid):
            file.write(f"{uniprotid} || {chemblid} || {get_bind(uniprotid)} || {get_activity(chemblid)}\n")
        else:
            continue
# 读取ID，并构造知识库

listt = []
with open('data_process/process.txt','r') as file:
    for i in file:
        list = i.split(' || ')
        uniprotid = list[0]
        chemblid = list[1]
        bind = list[2]
        activity = list[3]
        request = '……'
        S = promptt(pref_name=pref_name,Ki=activity,binding_site=bind,request=request,chemblid=chemblid)
        list = S.get_prompt()
        listt.append(list)

# 根据知识构造prompt

with open('data_process/prompt_list.json', 'w') as json_file:
    json.dump(listt, json_file, indent=4)

print(f"prompt的json已成功写入到 文件中。")

if __name__ == '__main__':
    with open('prompt_list.json', 'r') as jsonfile:
        dict_list = json.load(jsonfile)
    answer_list = []
    for i in dict_list:
        answer = get_answer([i[0],i[1],i[2],i[3]])
        print(answer)
        answer_list.append({"answer":answer,"chemblid":i[4]})

    with open('output/qwen_output.txt','w') as file:
        json.dump(answer_list, file, indent=4)

	pipe = pipeline("text-generation", model="/root/autodl-tmp/drugassist")
	with open('gpt4_shaixuan_jieguo.txt','r') as file:
    		Smiles = []
    		for i in file:
        	smiles = i.split(' || ')[0]
        	smiles = get_youhua(smiles,pipe)
        	Smiles.append(smiles+' || '+i.split(' || ')[1])

	with open('gpt4_shaixuan_jieguo.txt','w') as file:
    		for i in Smiles:
        	file.write(i+'\n')
#优化分子smiles并写回。
