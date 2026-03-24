# 自动处理 oil, protein, size 三个核心性状
for trait in oil protein size; do
    echo "--------------------------------------"
    echo "Processing trait: $trait"
    
    # 执行转换 (输入: .assoc.txt  输出: .ma  样本量: 378)
    python convert_gemma_to_gsmr.py \
        output/${trait}_analysis.assoc.txt \
        output/${trait}.ma \
        378
done

echo "--------------------------------------"
echo "Batch conversion complete. Check your .ma files."
ls -lh output/*.ma