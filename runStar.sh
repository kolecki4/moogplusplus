#! /bin/bash
ARGC=$#
#echo "███████████                                  _      _    "
#echo " ___    _     ______    ______    ______   _│ │_  _│ │_  "
#echo "│   \  / \   /  __  \  /  __  \  /  ____\ │_   _││_   _│ "
#echo "│  \ \/   \ │  /  \  ││  /  \  ││  /  ___   │_│    │_│   "
#echo "│  │\  /│  ││  │  │  ││  │  │  ││  │ │   │               "
#echo "│  │ \/ │  ││  \__/  ││  \__/  ││  \__/  │               "
#echo "│__│    │__│ \______/  \______/  \______/                "

echo " ╔═════════════════════════════════════════════════════════════════════════════╗"
echo " ║     *__________* ___    __  * ______    ______  * ______  *  _      _     * ║"
echo " ║     ___________ │   ╲  ╱  ╲  ╱  __  ╲ *╱  __  ╲  ╱  ____╲  _│ │_ *_│ │_     ║"
echo " ║    ____________ │  \ ╲╱    ││  ╱  ╲  ││  ╱  ╲  ││  ╱  ___ │_\`  _││_\`  _│    ║"
echo " ║ * _____________ │  │╲  ╱│  ││ │    │ ││ │    │ ││ │  │_  │  │_│    │_│      ║"
echo " ║  ______________ │  │ ╲╱ │  ││  ╲__╱  ││  ╲__╱  ││  ╲__╱  │    *     *       ║"
echo " ║ _______________ │__│ *  │__│ ╲______╱* ╲______╱  ╲______╱ *            *    ║"
echo " ║   *        *                 *                  *               *           ║"
echo " ╟─────────────────────────────────────────────────────────────────────────────╢"
echo " ║ ░░░░░░░▒▒▒▒▒▒▒▓▓▓▓▓▓▓████████ Beta Version 1.1 ███████▓▓▓▓▓▓▓▒▒▒▒▒▒▒░░░░░░░ ║"
echo " ╚═════════════════════════════════════════════════════════════════════════════╝"


if [ $ARGC -ne 3 ];
then
    echo "USAGE: ./runStar.sh [Star Name] [Working Directory] [Initial Metallicity]"
    echo ""    
    echo "  - Star Name: Must be resolvable by SIMBAD. 
    Used to query stellar photometry"
    echo ""
    echo "  - Working Directory: A full path to your star's working directory
    The final character of this path MUST BE A SLASH. This is 
    the folder in which your input data is located and where 
    output files will be dumped."
    echo ""    

    echo "  - Initial Metallicity: A \"starting guess\" for the code.
    If you have absolutely no idea, a value of 0 will work fine enough"
    
else
    starName=$1
    workDir=$2
    initMetal=$3

    i=0
    python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
    notConverged=1
    ((i++))
    until [[ $i -gt 0 && notConverged -eq 0 ]];
    do
        

        ./RunAbundanceOnGoodLines ${workDir}params.txt
        notConverged=$?
        python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
        ((i++))
    done
fi;
