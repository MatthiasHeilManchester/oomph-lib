# Setup directories
main_dir=BOTH_ADAPTATION_1pt2_again
if [ -e $main_dir ]; then
    echo "Directory " $main_dir "already exists -- please move it out of the way."
    exit
fi
mkdir $main_dir

# Mesh generators
#----------------
mesher_postfix_list="_gmsh _tetgen" #  _gmsh" # _tetgen"
code=""
for mesher_postfix in `echo $mesher_postfix_list`; do
    code=assess_tet_mesh_volume$mesher_postfix
    make $code
    cp $code $main_dir
    cd $main_dir
                
    run_dir=CASE`echo $mesher_postfix`
    mkdir $run_dir
    cp $code $run_dir
    cd $run_dir
    mkdir RESLT
    output_file=OUTPUT
    
    nref=10 # 4 # 5
    nunref=10 # 4 # 5

    target_volume=0.1
    if [ $mesher_postfix == "_gmsh" ]; then
        target_volume=1.0
    fi

    echo "Flags:  --suppress_inner_boundaries --desired_target_volume $target_volume --check_uniform_refinement $nref --check_uniform_unrefinement $nunref --specify_gmsh_target_volume_via_file --refinement_volume_decrease_factor 1.2"

    exit
    
    #./$code --suppress_inner_boundaries --desired_target_volume $target_volume --check_uniform_refinement $nref --check_uniform_unrefinement $nunref --specify_gmsh_target_volume_via_file --refinement_volume_decrease_factor 1.2  > $output_file


# --unrefinement_volume_increase_factor 3.0

#  --refinement_volume_decrease_factor 1.2  
    cd ..
    cd ..
done
