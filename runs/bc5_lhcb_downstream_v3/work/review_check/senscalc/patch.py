
import re
s=open('extract_nb3.py',encoding='utf-8').read()
old = """name=node[1]; args=node[2]
    a=[render(x) for x in args]"""
new = """name=node[1]; args=node[2]
    if name=='RowBox':
        if args and args[0][0]=='list':
            return ''.join(render(x) for x in args[0][1])
        return ''.join(render(x) for x in args)
    a=[render(x) for x in args]"""
assert old in s
s=s.replace(old,new)
open('extract_nb3.py','w',encoding='utf-8').write(s)
print('patched render')
