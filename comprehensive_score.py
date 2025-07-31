import math

def normalize_vina_score(vina_dock):
    """
    将Vina对接分数(负值)标准化为0-100评分系统
    分数越高代表结合能力越好
    
    参数:
        vina_dock: Vina对接分数(负值)，通常范围在-4到-12之间
        
    返回:
        float: 标准化分数(0-100)
    """
    # 将负值转换为正值 (x越大代表结合越好)
    x = -vina_dock
    
    # 基础S型曲线 (处理6-8.5区间)
    s_curve = 80 / (1 + math.exp(-4.0 * (x - 7.0)))
    
    # 优秀区增强函数 (x>8.5)
    excellence_boost = 0
    if x > 8.5:
        # 高斯增长函数 (先快后慢)
        excellence_boost = 20 * (1 - math.exp(-0.6 * (x - 8.5)**1.7))
    
    # 组合基础分数和优秀奖励
    base_score = s_curve + excellence_boost
    
    # 弱结合区惩罚 (x<6)
    if x < 6.0:
        # 指数衰减惩罚
        penalty = math.exp(4.0 * (x - 6.0))
        return max(0, base_score * penalty)
    
    return min(100, base_score)
def comprehensive_score(vina_dock, qed, sa):
    """
    计算分子综合质量评分
    :param vina_dock: 对接分数 (0-100)
    :param qed: 药物相似性 (0-1)
    :param sa: 合成可达性 (0-1)
    :return: 综合评分 (0-100)
    """
    # 限制vina_dock在[60,90]区间用于调整计算
    vina_dock = normalize_vina_score(vina_dock)
    clamped_v = max(60, min(vina_dock, 90))
    
    # 计算调整权重 (60/90处为0, 75处为1)
    weight = math.sin(math.pi * (clamped_v - 60) / 30)
    
    # 计算qed/sa的联合调整量
    qed_adj = 25 * (qed - 0.55)
    sa_adj = 25 * (sa - 0.55)
    total_adj = (qed_adj + sa_adj) * weight
    
    # 最终得分
    return vina_dock + total_adj