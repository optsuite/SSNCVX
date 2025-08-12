'''
Author: Hantao Nie (nht@pku.edu.cn)
Date: 2023-08-25 17:10:00
LastEditors: Hantao Nie (nht@pku.edu.cn)
LastEditTime: 2023-08-25 19:00:54
Description: 

Copyright (c) 2023, Hantao Nie, Peking University. 
'''
def begin_file():
    s = ""
    s += "<HTML>\n"
    s += "<HEAD>\n"
    s += "<TITLE>\n"
    s += "</TITLE>\n"
    s += "</HEAD>\n"
    s += "<BODY>\n"
    return s

def begin_table(table_name="table"):
    s = ""
    s += "<PRE>\n"
    s += "%20s" % table_name
    s += '\n'
    s += '*' * 250 + '\n'
    return s

def end_table():
    s = ""
    s += "\n"
    s += '*' * 250
    s += "\n"
    s += "</PRE>\n\n\n"
    return s


def end_file():
    s = ""
    s += "\n"
    s += '*' * 250
    s += "\n"
    s += "</PRE>\n\n\n"
    return s

    
def print(inp):
    lines = inp.split('\n')
    s = ""
    for line in lines:
        s += "<br>%s\n" % line
    return s

