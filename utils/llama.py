import openai
from openai import OpenAI

def chat(prompts, model="llama-3.1-405b", temperature=0.6):
    """
    多轮对话调用大模型并返回最后一轮的回复
    
    参数:
    prompts (list): 对话轮次的提示文本列表
    model (str): 要使用的模型名称，默认为"llama-3.1-405b"
    max_tokens (int): 生成的最大token数，默认为100
    temperature (float): 控制输出随机性的参数，默认为0.7
    
    返回:
    str: 最后一轮对话的模型回复
    """
    client = OpenAI(api_key = 'sk-meqtlZDPh2QLN7Je6e1b562cD99a4279B39985D417604f9b',base_url = 'https://api.expansion.chat/v1')
    messages = [
        {"role": "system", "content": "You are a helpful assistant."}
    ]
    
    for prompt in prompts:
        messages.append({"role": "user", "content": prompt})
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=temperature
        )
        messages.append({"role": "assistant", "content": response.choices[0].message.content.strip()})
    
    return messages[-1]["content"]
