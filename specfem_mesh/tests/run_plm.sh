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


echo "Running test_plm_function..."
./test_plm_function
echo "Running test_plm_values..."
./test_plm_values



if [ "$skip_pytest" = false ]; then
    pytest plm/test_plm.py
fi
if [ "$skip_cleanup" = false ]; then
    cd plm  
    bash cleanup.sh
fi