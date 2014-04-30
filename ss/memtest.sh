make clean; make; valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --show-reachable=yes --dsymutil=yes --xml=yes --xml-file="valtest.xml"  ./_ss.o;
python ../visualvalgrind/src/visualvalgrind.py build /Users/amabdol/Dropbox/Developer/GitHub/Banga_ScatterSearch/valtest.xml;

dot -Tsvg Leak_DefinitelyLost.dot > Leak_DefinitelyLost.svg;
dot -Tsvg Leak_IndirectlyLost.dot > Leak_IndirectlyLost.svg;
dot -Tsvg Leak_StillReachable.dot > Leak_StillReachable.svg;
dot -Tsvg UninitCondition.dot > UninitCondition.svg;
dot -Tsvg UninitValue.dot > UninitValue.svg;
dot -Tsvg InvalidFree.dot > InvalidFree.svg;
dot -Tsvg InvalidRead.dot > InvalidRead.svg;

# open -a safari Leak_DefinitelyLost.svg
# open -a safari Leak_IndirectlyLost.svg
open -a safari Leak_StillReachable.svg
# open -a safari UninitCondition.svg
# open -a safari UninitValue.svg
open -a safari InvalidFree.svg
open -a safari InvalidRead.svg