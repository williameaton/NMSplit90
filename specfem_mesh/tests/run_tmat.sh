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

echo "Running tmat_semi..."
./tmat_semi
echo "Running test_tmat..."
./test_tmat

if [ "$skip_pytest" = false ]; then
    pytest test_tmat_ic.py
fi
if [ "$skip_cleanup" = false ]; then
    cd Tmat/
    bash cleanup.sh
fi