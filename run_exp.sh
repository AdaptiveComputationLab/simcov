dir="experiments"
mkdir -p ${dir}

# j=(0.0007632 0.002049 0.004783 0.15 1)
j=(0.0007632)
# m=(0.007632 0.0717 0.763 1)
m=(0.763)
for valj in ${j[*]}; do
     for valm in ${m[*]}; do
        echo "${dir}/j-${valj}_m-${valm}"
        python create_config.py --virion-diffusion $valj --chemokine-diffusion $valm --sample-period 10080 --output "${dir}/j-${valj}_m-${valm}"  > "${dir}/j-${valj}_m-${valm}.config"
        sbatch run_p1.sh "${dir}/j-${valj}_m-${valm}.config"
    done
done