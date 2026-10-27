
import sys, re, io, os

def unescape(s):
    out=[]
    i=0
    while i < len(s):
        c=s[i]
        if c=='\\' and i+1<len(s):
            n=s[i+1]
            if n=='n': out.append('\n'); i+=2
            elif n=='t': out.append('\t'); i+=2
            elif n=='"': out.append('"'); i+=2
            elif n=='\\': out.append('\\'); i+=2
            elif n=='\n': out.append(''); i+=2   # line continuation
            else: out.append('\\'+n); i+=2
        else:
            out.append(c); i+=1
    return ''.join(out)

def strings(text):
    # returns list of (start, end, value) for double-quoted strings honoring escapes
    res=[]; i=0; n=len(text)
    while i<n:
        if text[i]=='"':
            j=i+1; buf=[]
            while j<n:
                if text[j]=='\\':
                    buf.append(text[j]); 
                    if j+1<n: buf.append(text[j+1])
                    j+=2
                elif text[j]=='"':
                    break
                else:
                    buf.append(text[j]); j+=1
            res.append((i,j,unescape(''.join(buf))))
            i=j+1
        else:
            i+=1
    return res

def main(path):
    text=open(path,encoding='utf-8',errors='replace').read()
    # split into cells
    cells=text.split('Cell[')
    print('### FILE:', path)
    for k,cell in enumerate(cells[1:],1):
        ss=strings(cell)
        if not ss: continue
        joined=''.join(v for _,_,v in ss)
        if not joined.strip(): continue
        print('\n----- cell %d -----'%k)
        print(joined)

main(sys.argv[1])
