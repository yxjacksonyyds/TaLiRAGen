from test import test
import multiprocessing
from multiprocessing import Process, Queue
from contextlib import contextmanager
@contextmanager
def timeout_managed(func, timeout, *args, **kwargs):
    queue = Queue()
    def wrapper():
        try:
            result = func(*args, **kwargs)
            queue.put(result)
        except Exception as e:
            queue.put(e)
    
    p = Process(target=wrapper)
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        p.join()
        yield None  # 超时返回 None
    else:
        result = queue.get()
        if isinstance(result, Exception):
            raise result
        else:
            yield result

from rdkit import Chem

def is_valid_smiles(smiles):
    """
    判断给定的SMILES字符串是否合理。

    :param smiles: str, SMILES字符串
    :return: bool, 如果SMILES字符串合理则返回True，否则返回False
    """
    try:
        # 尝试从SMILES字符串创建一个分子对象
        molecule = Chem.MolFromSmiles(smiles)
        # 如果分子对象创建成功并且不是None，则SMILES字符串是合理的
        return molecule is not None
    except:
        # 如果在创建分子对象过程中出现异常，则SMILES字符串不合理
        return False

def get_score(smiles,pref_name):
    t = test()
    qed,sa,affinity = t.cauculate(smiles,pref_name)
    return qed,sa,affinity

def avg(float_list):
    total_sum = sum(float_list)
    
    # 计算列表中值的数量
    count = len(float_list)
    
    # 计算平均数
    average = total_sum / count
    return average

def deal_file(path):
    Qed = []
    Affinity = []
    SA = []
    time = 0
    valid = 0
    with open(path,'r') as file:
        for line in file:
            try:
                print("1")
                sentence = line.split(' || ')
                smiles = sentence[0]
                pref_name = sentence[1][:-1]
                # 使用带超时的上下文管理器执行get_score
                if is_valid_smiles(smiles):
                    with timeout_managed(get_score, 120, smiles, pref_name) as result:
                        if result is None:
                            print(f"Timeout for {smiles}, skipping...")
                            continue
                        qed, sa, affinity = result  # 正确解包结果
                        if affinity < -5:
                            Qed.append(qed)
                            SA.append(sa)
                            Affinity.append(affinity)
                            valid = valid +1
                            print(time,avg(Qed),avg(Affinity),avg(SA),valid)

                else:
                    print('invalid smiles')
            except Exception:
                continue
            
            time = time + 1
    return avg(Qed),avg(Affinity),avg(SA)

if __name__ == '__main__':
    print(deal_file('/root/autodl-tmp/gpt4_shaixuan_jieguo.txt'))
        