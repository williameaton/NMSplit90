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


echo "Running test_xlm..."
./test_xlm


if [ "$skip_pytest" = false ]; then
    pytest xlm/test_xlm.py
fi
if [ "$skip_cleanup" = false ]; then
    cd xlm 
    bash cleanup.sh
fi