#!/usr/bin/env python3
import re
import shutil
from pathlib import Path

cp_pattern = r'cp\s+(-[RHLfiPnaclpSsvXx]*)?\s*([\w\.\-\/]+)\s+([\w\.\-\/]+)'

# 定义cp命令的正则表达式模式
def execute_cp_command(command):
    match = re.match(cp_pattern, command)
    if match:
        options, source, destination = match.groups()
        options = options if options else "No options"
        print(f"Options: {options}, Source: {source}, Destination: {destination}")
        
        # 这里只是一个示意，实际复制文件时需要根据选项进行相应操作
        try:
            # shutil.copy(source, destination)
            print(f"File copied from {source} to {destination}")
        except Exception as e:
            print(f"Error copying file: {e}")
    else:
        print("Invalid cp command syntax.")

# 示例使用
# execute_cp_command("cp source_file.txt target_file.txt")
# execute_cp_command("cp -R source_directory target_directory")
import re

def resolve_placeholders(data, text, separator=':'):
    # 定义一个正则表达式来匹配形式为 ${...} 的模式
    pattern = re.compile(r'\$\{([^}]+)\}')

    # 定义一个递归函数来根据键路径查找值
    def find_value(d, keys):
        if keys and keys[0] in d:
            return find_value(d[keys[0]], keys[1:]) if len(keys) > 1 else d[keys[0]]
        return None

    # 查找所有匹配的模式
    for match in pattern.finditer(text):
        # 获取键路径，如 'option:Sub'
        key_path = match.group(1)
        keys = key_path.split(separator)
        # 根据键路径查找值
        value = find_value(data, keys)
        # 如果找到了值，替换原文中的模式
        if value is not None:
            text = text.replace(match.group(0), value)

    return text

# 示例数据和文本
data = {"option": {"Sub": "some infos", "Class": {"type": "MORE USEFUL TAGGING"}}}
text = "This is a test string with ${option:Sub} to be replaced."

# 使用函数进行替换
result_text = resolve_placeholders(data, text)
print(result_text)

text = "This is a test string with ${option:Class:type} to be replaced."

# 使用函数进行替换
result_text = resolve_placeholders(data, text)
print(result_text)