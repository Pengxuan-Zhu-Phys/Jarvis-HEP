import sympy as sp

def evaluate_expression(expr_str, variables):
    """
    使用 sympy 对表达式进行求值。

    参数:
        expr_str (str): 表达式字符串，例如 "(2.0 * X < Y)"。
        variables (dict): 变量名和对应值的字典，例如 {"X": 5.0, "Y": 10}。

    返回:
        bool: 表达式是否为真。
    """
    # 将变量转换为 sympy 符号
    symbols = {var: sp.symbols(var) for var in variables}
    
    # 将字符串解析为 sympy 表达式
    expr = sp.sympify(expr_str, locals=symbols)
    
    # 用给定变量的值进行求值
    result = expr.subs(variables)
    
    # 返回布尔值结果
    return bool(result)

# 示例用法
expression = "(2.0 * X < Y)"
variables = {"X": 5.0, "Y": 10}
result = evaluate_expression(expression, variables)
print(result)  # 输出: True