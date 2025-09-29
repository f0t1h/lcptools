
echo "#ifndef _LCPTOOLS_HO_H_" > lcptools_ho.h
echo "#define _LCPTOOLS_HO_H_" >> lcptools_ho.h
for hfile in lps.h encoding.h core.h; do
    cat $hfile | grep -v "#include \""  >> lcptools_ho.h
done
echo "#ifdef LCPTOOLS_IMPL"  >> lcptools_ho.h
for cfile in lps.c encoding.c core.c; do
    cat $cfile | grep -v "#include \""  >> lcptools_ho.h
done
echo "#endif" >> lcptools_ho.h
echo "#endif" >> lcptools_ho.h
