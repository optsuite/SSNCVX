#!/bin/bash
#SBATCH -J R1TA      # 任务名字是 matlab-test
#SBATCH -o R1TA.out     # 输出文件是 matlab-test.out
#SBATCH -p cpu                 # 提交到 cpu 分区
#SBATCH -w c04             # 提交到 cpu 分区
#SBATCH -N 1                   # 申请 1 节点
#SBATCH --cpus-per-task=4      # 申请 4 cpu core
#SBATCH --cpus-per-cpu=50000 
#SBATCH -t 7-00:00:00               # 最长运行时间 30 分钟


module add matlab              # 添加 MATLAB 模块
# 使用 MATLAB 运行当前文件夹中的 matlab-test.m 文件
matlab -batch "Test_R1TA"

# 对 R2019a 及之后版本的 MATLAB 还有更简洁的用法
# 取消注释可以看到效果
#matlab -batch "matlab_test"