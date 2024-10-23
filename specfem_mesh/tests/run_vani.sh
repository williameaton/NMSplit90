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


echo "Running test_vani_diag..."
./test_vani_diag
echo "Running test_vani_fast..."
./test_vani_fast
echo "Running test_tromp1995..."
./test_tromp1995


if [ "$skip_pytest" = false ]; then
    pytest Vani/test_vani.py
fi
if [ "$skip_cleanup" = false ]; then
    cd Vani/
    bash cleanup.sh
fi