#!/bin/bash
#SBATCH -J ssncp-socp-matlab         # 任务名字是 matlab-test
#SBATCH -p cpu                       # 提交到 cpu 分区
#SBATCH -N 1                         # 申请 1 个节点
#SBATCH --cpus-per-task=64           # 申请 4 个 CPU 核心
#SBATCH -t 24:00:00                  # 最长运行时间 6 小时

# 加载 Matlab 模块
module unload matlab
module load matlab/r2023a
export LD_PRELOAD=/bicmr/soft/MATLAB/R2023a/bin/glnxa64/libmwmetis.so

# 运行 Matlab 脚本
matlab -nodesktop -nosplash -nodisplay -r "Test_mittelmann"