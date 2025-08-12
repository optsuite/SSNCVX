'''
Author: Hantao Nie (nht@pku.edu.cn)
Date: 2023-08-25 13:52:25
LastEditors: Hantao Nie (nht@pku.edu.cn)
LastEditTime: 2023-08-25 21:27:56
Description: 

Copyright (c) 2023, Hantao Nie, Peking University. 
'''
import numpy as np
import pandas as pd
import os
import math
import scipy.io

import html_utils as html


string_format = "%10s"
long_string_format = "%20s"
long_long_string_format = "%30s"
float_format = "%10.2f"
percent_format = "%9d%%"
scientific_format = "%10.1e"
int_format = "%10d"


def format_print(value, is_ratio=False, 
            format="html",
            prefer_scitific=False,
            prefer_int=True) -> str:
    if not isinstance(value, float):
        return string_format % value
    if value == -1:  #表示空位，用在iter等指标的geomean
        return string_format % ""
    if value == -2:  #表示读取fail
        return string_format % "--"
    if is_ratio:
        out = percent_format % (value * 100)
        if format == "latex":
            out = out.replace("%", "\\%")
        return out
    if prefer_scitific or abs(value) >= 1e4:
        return scientific_format % value
    if abs(value - round(value)) < 1e-16: #detect value as integer
        return int_format % value
    if abs(value) <= 1e-3:
        return scientific_format % value

    if prefer_int and abs(value) > 10:
        return int_format % round(value)
    return float_format % value

def read_txt(file_name, instance, method, info_list):
    out = {}
    try:
        f = open(file_name, 'r')
        result = f.read().split()
        f.close()
        for i in range(len(result)):
            x = result[i]
            try:
                x = float(x) 
                if math.isnan(x) or x >= 1e10:
                    x = -1
            except:
                x = str(x)
            out[method + "." + info_list[i]] = x

        if len(result) < len(info_list):
            for i in range(len(result), len(info_list)):
                out[method + "." + info_list[i]] = -1
            print("results is incomplete", instance, method)
        elif len(result) > len(info_list):
            print("results is over complete", instance, method)
    except:
        print("read result failed ", instance, method)
        for info in info_list:
            out[method + "." + info] = -2
    return out


def read_mat(file_name, instance, method, info_list):
    print(file_name)
    out = {}
    try:
        mat = scipy.io.loadmat(file_name)
        for info in info_list:
            x = mat.get(info)
            if x is None:
                out[method + "." + info] = -1
                print("results is incomplete", instance, method)
            else:
                try:
                    x = float(x)
                    if math.isnan(x) or x >= 1e10:
                        x = -1
                except:
                    x = str(x)
                out[method + "." + info] = x
    except:
        print("read result failed ", instance, method)
        for info in info_list:
            out[method + "." + info] = -2
    return out

class MyTable:
    def __init__(self, name="table", instance_list=None, method_list=None, info_list=None):
        self.name = name

        self.instance_list = instance_list
        self.method_list = method_list
        self.info_list = info_list

        
        # intialize the table
        self.data = pd.DataFrame(columns=[method + "." + info for method in method_list for info in info_list])
        
        for instance in instance_list:
            self.data.loc[instance] = [-1 for _ in range(len(self.data.columns))]


    def load_data(self, data_dir, read_file_func = None, filename_style = None):
        if read_file_func is None:
            read_file_func = read_txt
        if filename_style is None:
            filename_style = lambda instance, method: os.path.join(instance, method + ".txt")
        
        for instance in self.instance_list:
            for method in self.method_list:
                file_name = os.path.join(data_dir, filename_style(instance, method))
                result = read_file_func(file_name, instance, method, self.info_list)
                self.data.loc[instance, result.keys()] = list(result.values())
        


    def compute_acc(self):
        for method in self.method_list:
            self.data[method + ".acc"] = self.data[[method + ".gap", method + ".pinf", method + ".dinf"]].max(axis=1) 

    def compute_fail_penalty(self, info_list = ["time"], option = 3600):
        assert option in ["2maxrow", "3maxrow"] or isinstance(option, int) or isinstance(option, float)
        for instance in self.instance_list:
            for info in info_list:
                if option == "2maxrow":
                    cols = [method + "." + info for method in self.method_list]
                    self.data.loc[instance, "fail_penalty" + "." + info] = self.data.loc[instance, cols].max() * 2
                elif option == "3maxrow":
                    cols = [method + "." + info for method in self.method_list]
                    self.data.loc[instance, "fail_penalty" + "." + info] = self.data.loc[instance, cols].max() * 3
                else:
                    self.data.loc[instance, "fail_penalty" + "." + info] = option

    
    def compute_mean(self, option = "geomean", info_list = ["time"], shift=10):
        assert option in ["geomean", "mean"]
        self.data.loc[option] = [-1 for _ in range(len(self.data.columns))]
        for method in self.method_list:
            for info in info_list:
                if not "fail_penalty" + "." + info in self.data.columns:
                    self.compute_fail_penalty(info_list = [info])
                col = method + "." + info
                filtered_data = self.data.loc[(self.data[col] >= 0) & (self.data.index.isin(self.instance_list)), col]
                fail_penalty_data = self.data.loc[(self.data[col] < 0) & (self.data.index.isin(self.instance_list)), "fail_penalty" + "." + info]                
                filtered_data = pd.concat([filtered_data, fail_penalty_data])
                if filtered_data.empty:
                    mean_value = -1
                else:
                    if option == "mean":
                        mean_value = filtered_data.mean()
                    elif option == "geomean":
                        mean_value = math.exp(np.log(filtered_data + shift).mean()) - shift

                self.data.loc[option, col] = mean_value


    def set_fail(self, instance, method):
        # set all the info to -2
        for info in self.info_list:
            self.data.loc[instance][method + "." + info] = -2

    def check_acc(self, tol):
        for instance in self.instance_list:
            for method in self.method_list:
                if method + ".acc" not in self.data.columns:
                    self.compute_acc()
                if self.data.loc[instance][method + ".acc"] < 0 or self.data.loc[instance][method + ".acc"] > tol:
                    print("acc check failed", instance, method)
                    self.set_fail(instance, method)

    def check_status(self):
          for instance in self.instance_list:
            for method in self.method_list:
                if self.data.loc[instance][method + ".status"] != 's':
                    print("status check failed", instance, method)
                    self.set_fail(instance, method)

    def export(self, outfile_path = None, rows = None, method_list = None, info_list = None, format="html"):
        if rows is None:
            rows = self.data.index
        if method_list is None:
            method_list = self.method_list
        if info_list is None:
            info_list = self.info_list
        assert format in ["html", "latex"]

        if format == "html":
            endl = "\n"
            sep = ""
            row_name_format = long_long_string_format
        else:
            #endl="\\\\\hline\n"
            endl="\\\\\n"
            sep ="&"
            row_name_format = "%s"
            
        s = ""
        if format == "html":
            s += html.begin_file()
            s += html.begin_table(table_name=self.name)

        # generate the header
        if format == "html":
            s += row_name_format % ""
            for method in method_list:
                s += sep + string_format % method
                for _ in range(len(self.info_list) - 1):
                    s += sep + string_format % ""
            s += endl

            s += row_name_format %  "problem"
            for method in self.method_list:
                for info in self.info_list:
                    s += sep + string_format % info
            s += endl
            s += '*' * 250
            s += endl

        # generate the body
        for i in range(len(rows)):
            row = rows[i]
            if format == "latex" and row in ["mean", "geomean"]:
                s += "\\hline\n"
            s += row_name_format % (row.replace("_", "\\_") if format == "latex" else row)
            prefer_int = not (row ==  "geomean" )
            for method in method_list:
                for info in info_list:
                    prefer_scitific = info == "acc"
                    x = self.data.loc[row][method + "." + info]
                    s += sep + format_print(x, is_ratio=False, format=format, prefer_scitific=prefer_scitific, prefer_int=prefer_int)             
            if i < len(rows)-1:
                s += endl         
        if format == "html":
            s += html.end_table()
        if outfile_path is not None:
            f = open(outfile_path, 'w')
            f.write(s)
            f.close()
        return s