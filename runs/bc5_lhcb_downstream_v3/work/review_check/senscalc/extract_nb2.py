
import sys, re

def unescape(s):
    out=[]; i=0
    while i < len(s):
        c=s[i]
        if c=='\\' and i+1<len(s):
            n=s[i+1]
            if n=='n': out.append('\n'); i+=2
            elif n=='t': out.append('\t'); i+=2
            elif n=='"': out.append('"'); i+=2
            elif n=='\\': out.append('\\'); i+=2
            elif n=='\n': i+=2
            else: out.append('\\'+n); i+=2
        else:
            out.append(c); i+=1
    return ''.join(out)

def tokenize(text):
    toks=[]; i=0; n=len(text)
    while i<n:
        c=text[i]
        if c in ' \t\r\n': i+=1; continue
        if c=='"':
            j=i+1; buf=[]
            while j<n:
                if text[j]=='\\':
                    buf.append(text[j]);
                    if j+1<n: buf.append(text[j+1])
                    j+=2
                elif text[j]=='"': break
                else: buf.append(text[j]); j+=1
            toks.append(('str', unescape(''.join(buf)))); i=j+1; continue
        if c in '[]{},':
            toks.append((c,c)); i+=1; continue
        # identifier: letters, digits, $, backtick, underscore, ^ (caret used in symbols?), *? 
        m=re.match(r'[A-Za-z0-9$\x60.^\\]+', text[i:])
        if m:
            toks.append(('id', m.group(0))); i+=len(m.group(0)); continue
        toks.append(('other', c)); i+=1
    return toks

class P:
    def __init__(self,toks): self.t=toks; self.i=0
    def peek(self): return self.t[self.i] if self.i<len(self.t) else (None,None)
    def next(self):
        v=self.t[self.i]; self.i+=1; return v
    def parse(self):
        typ,val=self.peek()
        if typ=='str': self.next(); return ('str',val)
        if typ in '{}':
            return self.parse_list()
        if typ=='id':
            self.next()
            if self.peek()[0]=='[':
                self.next()
                args=[]
                if self.peek()[0]==']':
                    self.next()
                else:
                    while True:
                        args.append(self.parse())
                        t,v=self.next()
                        if t==']': break
                        if t!=',': 
                            # unexpected, bail
                            pass
                return ('call',val,args)
            return ('id',val)
        # other
        self.next(); return ('other',val)
    def parse_list(self):
        self.next() # {
        items=[]
        if self.peek()[0]=='}':
            self.next(); return ('list',items)
        while True:
            items.append(self.parse())
            t,v=self.next()
            if t=='}': break
        return ('list',items)

def render(node):
    k=node[0]
    if k=='str': return node[1]
    if k=='id': return node[1]
    if k=='other': return node[1]
    if k=='list': return '{'+', '.join(render(x) for x in node[1])+'}'
    if k=='call':
        name=node[1]; args=node[2]
        a=[render(x) for x in args]
        if name=='RowBox':
            return ''.join(a[0].strip('{}').split(', ')) if a else ''
        if name=='SuperscriptBox': return '%s^%s'%(a[0],a[1]) if len(a)>1 else a[0]
        if name=='FractionBox': return '(%s)/(%s)'%(a[0],a[1]) if len(a)>1 else a[0]
        if name=='SubscriptBox': return '%s_%s'%(a[0],a[1]) if len(a)>1 else a[0]
        if name=='SqrtBox': return 'Sqrt[%s]'%(a[0],)
        if name=='RadicalBox': return '(%s)^(1/%s)'%(a[0],a[1])
        if name=='StyleBox': return a[0]
        if name=='TagBox': return a[0]
        if name=='FormBox': return a[0]
        if name=='InterpretationBox': return a[0]
        if name=='TooltipBox': return a[0]
        if name=='ButtonBox': return a[0]
        if name=='UnderscriptBox': return 'Underscript[%s,%s]'%(a[0],a[1])
        if name=='OverscriptBox': return '%s_(%s)'%(a[0],a[1])
        if name=='GridBox': return '<GRID>'
        return '%s[%s]'%(name,', '.join(a))
    return str(node)

def extract(path):
    text=open(path,encoding='utf-8',errors='replace').read()
    # find cells by scanning for 'Cell[' at top level-ish
    idxs=[m.start() for m in re.finditer(r'(?m)^Cell\[', text)]
    print('### FILE:', path)
    for n,st in enumerate(idxs):
        # find matching bracket for the Cell[
        i=st+4  # position of '['
        depth=0; j=i
        instr=False
        while j < len(text):
            c=text[j]
            if instr:
                if c=='\\': j+=2; continue
                if c=='"': instr=False
                j+=1; continue
            if c=='"': instr=True; j+=1; continue
            if c in '[{': depth+=1
            elif c in ']}':
                depth-=1
                if depth==0: break
            j+=1
        celltext=text[st:j+1]
        # find BoxData[ ... ] inside
        bi=celltext.find('BoxData[')
        if bi<0: continue
        # find matching bracket
        i2=bi+len('BoxData')  # at '['
        depth=0; k=i2; instr=False
        while k<len(celltext):
            c=celltext[k]
            if instr:
                if c=='\\': k+=2; continue
                if c=='"': instr=False
                k+=1; continue
            if c=='"': instr=True; k+=1; continue
            if c in '[{': depth+=1
            elif c in ']}':
                depth-=1
                if depth==0: break
            k+=1
        boxtxt=celltext[i2+1:k]
        try:
            toks=tokenize(boxtxt)
            p=P(toks); node=p.parse()
            out=render(node)
        except Exception as e:
            out='<PARSE ERROR %s: %s>'%(e, boxtxt[:200])
        # cell label line: find 'Cell["' label
        lab=re.search(r'Cell\["((?:[^"\\\\]|\\\\.)*)"', celltext)
        print('\n----- cell %d  label=%s -----'%(n+1, lab.group(1) if lab else ''))
        print(out)

extract(sys.argv[1])
