
# Smooth only
bash /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/WheatComp.sh HMM_smoother \
    -i /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/data \
    --folder_lis /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/folders.txt \
    -o /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/out\
    --processes 21

# Train & Smooth
bash /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/WheatComp.sh HMM_smoother \
    -i /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/data \
    --folder_lis /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/folders.txt \
    -o /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/out \
    --processes 21 \
    --train \
    --niter 30