
import sys, re

CHAR = {
 'IndentingNewLine':'\n','Rule':'->','RuleDelayed':':>','Equal':'==','LessEqual':'<=',
 'GreaterEqual':'>=','NotEqual':'!=','ThinSpace':'','VeryThinSpace':'','NonBreakingSpace':'',
 'Times':'*','Divide':'/','Minus':'-','CapitalUpsilon':'Upsilon','Pi':'Pi','Theta':'theta',
 'Epsilon':'epsilon','Gamma':'gamma','Alpha':'alpha','Beta':'beta','Mu':'mu','Tau':'tau',
 'LessLess':'','GreaterGreater':'','LeftAngleBracket':'<','RightAngleBracket':'>',
 'InvisibleSpace':'','InvisibleComma':'','NoBreak':'','AlignmentMarker':'',
 'LeftDoubleBracket':'[[','RightDoubleBracket':']]','Transpose':'','InvisibleApplication':'',
 'VeryThinSpace':''
}

def unescape(s):
    out=[]; i=0
    while i < len(s):
        c=s[i]
        if c=='\\' and i+1<len(s):
            n=s[i+1]
            if n=='<': out.append(''); i+=2
            elif n=='>': out.append(''); i+=2
            elif n=='n': out.append('\n'); i+=2
            elif n=='t': out.append('\t'); i+=2
            elif n=='"': out.append('"'); i+=2
            elif n=='\\': out.append('\\'); i+=2
            elif n=='\n': i+=2
            elif n=='[':
                m=re.match(r'\\\[([A-Za-z]+)\]', s[i:])
                if m:
                    nm=m.group(1)
                    out.append(CHAR.get(nm,'\\['+nm+']'))
                    i+=len(m.group(0))
                else:
                    out.append('\\'); i+=1
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
                    buf.append(text[j])
                    if j+1<n: buf.append(text[j+1])
                    j+=2
                elif text[j]=='"': break
                else: buf.append(text[j]); j+=1
            toks.append(('str', unescape(''.join(buf)))); i=j+1; continue
        if c in '[]{},':
            toks.append((c,c)); i+=1; continue
        m=re.match(r'[A-Za-z0-9$\x60.^]+', text[i:])
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
        if typ in '{}': return self.parse_list()
        if typ=='id':
            self.next()
            if self.peek()[0]=='[':
                self.next(); args=[]
                if self.peek()[0]==']': self.next()
                else:
                    while True:
                        args.append(self.parse())
                        t,v=self.next()
                        if t==']': break
                return ('call',val,args)
            return ('id',val)
        self.next(); return ('other',val)
    def parse_list(self):
        self.next()
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
    if k in ('str','id'): return node[1]
    if k=='other': return node[1]
    if k=='list': return '{'+', '.join(render(x) for x in node[1])+'}'
    name=node[1]; args=node[2]
    if name=='RowBox':
        if args and args[0][0]=='list':
            return ''.join(render(x) for x in args[0][1])
        return ''.join(render(x) for x in args)
    a=[render(x) for x in args]
    if name=='RowBox':
        if args and args[0][0]=='list':
            return ''.join(render(x) for x in args[0][1])
        return ''.join(a)
    if name=='SuperscriptBox': return '%s^%s'%(a[0],a[1]) if len(a)>1 else a[0]
    if name=='FractionBox': return '(%s)/(%s)'%(a[0],a[1]) if len(a)>1 else a[0]
    if name=='SubscriptBox': return '%s_%s'%(a[0],a[1]) if len(a)>1 else a[0]
    if name=='SqrtBox': return 'Sqrt[%s]'%(a[0],)
    if name=='RadicalBox': return '(%s)^(1/%s)'%(a[0],a[1])
    if name in ('StyleBox','TagBox','FormBox','InterpretationBox','TooltipBox','ButtonBox','ItemBox','PaneBox','AdjustmentBox','BoxFormBox','PanelBox','FrameBox'):
        return a[0] if a else ''
    if name=='UnderscriptBox': return 'Underscript[%s,%s]'%(a[0],a[1])
    if name=='OverscriptBox': return '%s_(%s)'%(a[0],a[1])
    if name=='GridBox': return '<GRID>'
    if name=='BoxData': return a[0]
    return '%s[%s]'%(name,', '.join(a))

def find_matching(text, i):
    # i points at opening bracket
    depth=0; j=i; instr=False
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
            if depth==0: return j
        j+=1
    return len(text)-1

def extract(path):
    text=open(path,encoding='utf-8',errors='replace').read()
    idxs=[m.start() for m in re.finditer(r'(?m)^Cell\[', text)]
    print('### FILE:', path)
    n=0
    for st in idxs:
        i=st+4
        j=find_matching(text,i)
        celltext=text[st:j+1]
        bi=celltext.find('BoxData[')
        lab=re.search(r'Cell\["((?:[^"\\\\]|\\\\.)*)"', celltext)
        label=unescape(lab.group(1)) if lab else ''
        if bi<0:
            continue
        i2=bi+len('BoxData')
        k=find_matching(celltext,i2)
        boxtxt=celltext[i2+1:k]
        n+=1
        try:
            toks=tokenize(boxtxt)
            p=P(toks); node=p.parse()
            out=render(node)
        except Exception as e:
            out='<PARSE ERROR %s: %s>'%(e, boxtxt[:300])
        print('\n----- cell %d  label=%s -----'%(n, label))
        print(out)

extract(sys.argv[1])
