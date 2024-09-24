# Check for flags
skip_cleanup=false
skip_pytest=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --nc) skip_cleanup=true ;;
        --np) skip_pytest=true ;;
    esac
    shift
done


echo "Running test_wrot..."
./test_wrot
echo "Running test_wstored..."
./test_wstored
echo "Running rotation_semi..."
./rotation_semi



if [ "$skip_pytest" = false ]; then
    pytest rot_mat/test_rotation_matrix.py
fi
if [ "$skip_cleanup" = false ]; then
    cd rot_mat
    bash cleanup.sh
fi