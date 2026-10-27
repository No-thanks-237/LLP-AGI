
import re
s=open('extract_nb3.py',encoding='utf-8').read()
old="out.append(CHAR.get(nm,'['+nm+']'))"
new="out.append(CHAR.get(nm,'\\\\['+nm+']'))"
assert old in s, 'not found'
s=s.replace(old,new)
open('extract_nb3.py','w',encoding='utf-8').write(s)
print('patched named-char default')
