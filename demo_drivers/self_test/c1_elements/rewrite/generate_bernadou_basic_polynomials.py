#!/usr/bin/env python3
"""Generate loop-free Bernadou basic-polynomial evaluators from C++ matrix entries.

Usage:
  python3 generate_bernadou.py SOURCE.cc 3 7 > bernadou_basic_3.h
  python3 generate_bernadou.py SOURCE.cc 5 9 > bernadou_basic_5.h

SOURCE.cc may contain the entire monomial_to_basic_matrix() definition. Only
assignments of the form A(i,j)=...; or b2l(i,j)=...; are read.
"""
import argparse, ast, math, re
from fractions import Fraction

ASSIGN = re.compile(r'\b(?:A|b2l)\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*(.*?);')

# Restricted exact evaluator: rational arithmetic plus sqrt(2), represented as a+b*sqrt(2).
class Q2:
    def __init__(self, a=0, b=0): self.a, self.b = Fraction(a), Fraction(b)
    def __add__(self,o): o=q2(o); return Q2(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self): return Q2(-self.a,-self.b)
    def __sub__(self,o): return self+(-q2(o))
    def __rsub__(self,o): return q2(o)-self
    def __mul__(self,o):
        o=q2(o); return Q2(self.a*o.a+2*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=q2(o); d=o.a*o.a-2*o.b*o.b
        return Q2((self.a*o.a-2*self.b*o.b)/d,(self.b*o.a-self.a*o.b)/d)
    def __rtruediv__(self,o): return q2(o)/self
    def iszero(self): return self.a==0 and self.b==0

def q2(x): return x if isinstance(x,Q2) else Q2(x)

def eval_expr(text):
    text=text.replace('std::sqrt(2.)','SQRT2').replace('std::sqrt(2.0)','SQRT2').replace('sqrt(2)','SQRT2').replace('sqrt(2.)','SQRT2')
    node=ast.parse(text,mode='eval').body
    def ev(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,(int,float)):
            return Q2(Fraction(str(n.value)))
        if isinstance(n,ast.Name) and n.id=='SQRT2': return Q2(0,1)
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub): return -ev(n.operand)
        if isinstance(n,ast.BinOp):
            x,y=ev(n.left),ev(n.right)
            if isinstance(n.op,ast.Add): return x+y
            if isinstance(n.op,ast.Sub): return x-y
            if isinstance(n.op,ast.Mult): return x*y
            if isinstance(n.op,ast.Div): return x/y
        raise ValueError(f'Unsupported coefficient expression: {ast.dump(n)}')
    return ev(node)

def monomials(n):
    return [(l-m,m) for l in range(n,-1,-1) for m in range(l+1)]

def frac_cpp(x):
    if x.denominator==1: return f'{x.numerator}.0'
    return f'({x.numerator}.0 / {x.denominator}.0)'

def coeff_cpp(c):
    if c.b==0: return frac_cpp(c.a)
    if c.a==0:
        if c.b==1: return 'sqrt_two'
        if c.b==-1: return '-sqrt_two'
        return f'{frac_cpp(c.b)} * sqrt_two'
    return f'({frac_cpp(c.a)} + {frac_cpp(c.b)} * sqrt_two)'

def varpow(v,n): return '1.0' if n==0 else v if n==1 else f'{v}{n}'
def mono(a,b):
    z=[]
    if a: z.append(varpow('s',a))
    if b: z.append(varpow('t',b))
    return ' * '.join(z) if z else '1.0'

def signed_terms(terms):
    if not terms: return ['0.0']
    return [f'{coeff_cpp(c)} * {mono(a,b)}' for c,a,b in terms]

def emit_assignment(lhs,terms,indent='  '):
    xs=signed_terms(terms)
    if len(xs)==1: return [f'{indent}{lhs} = {xs[0]};']
    out=[f'{indent}{lhs} =']
    for k,x in enumerate(xs): out.append(f'{indent}  {x}{";" if k==len(xs)-1 else " +"}')
    return out

def differentiated(row,axis0,axis1):
    out=[]
    for c,a,b in row:
        if a<axis0 or b<axis1: continue
        f=1
        for k in range(axis0): f*=a-k
        for k in range(axis1): f*=b-k
        out.append((c*f,a-axis0,b-axis1))
    return out

def function(bo,n,name,argtype,outname,derivs,rows):
    out=['template<>',f'void BernadouElementBasis<{bo}>::{name}(',f'  const Vector<double>& s_basic, {argtype}& {outname}) const','{','  const double s = s_basic[0];','  const double t = s_basic[1];']
    maxdeg=max(n-dx-dy for dx,dy in derivs)
    for k in range(2,maxdeg+1): out.append(f'  const double s{k} = s{k-1} * s;' if k>2 else '  const double s2 = s * s;')
    for k in range(2,maxdeg+1): out.append(f'  const double t{k} = t{k-1} * t;' if k>2 else '  const double t2 = t * t;')
    if any(c.b for r in rows for c,_,_ in r): out.append('  static const double sqrt_two = std::sqrt(2.0);')
    out.append('')
    for i,row in enumerate(rows):
        for col,(dx,dy) in enumerate(derivs):
            lhs=f'{outname}[{i}]' if len(derivs)==1 else f'{outname}({i}, {col})'
            out.extend(emit_assignment(lhs,differentiated(row,dx,dy)))
    out+=['}','']
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('source'); ap.add_argument('boundary_order',type=int); ap.add_argument('basis_order',type=int); args=ap.parse_args()
    text=open(args.source).read(); entries={}
    for i,j,e in ASSIGN.findall(text): entries[int(i),int(j)]=eval_expr(e)
    nmono=(args.basis_order+1)*(args.basis_order+2)//2
    if not entries: raise SystemExit('No matrix assignments found')
    if max(max(i,j) for i,j in entries)>=nmono: raise SystemExit('Matrix index exceeds expected size')
    powers=monomials(args.basis_order)
    rows=[]
    for i in range(nmono):
        rows.append([(entries[i,j],*powers[j]) for j in range(nmono) if (i,j) in entries and not entries[i,j].iszero()])
    out=['// Automatically generated from monomial_to_basic_matrix().','// Matrix coefficients were parsed as exact elements of Q(sqrt(2)).','']
    out += function(args.boundary_order,args.basis_order,'full_basic_polynomials','Shape','phi',[(0,0)],rows)
    out += function(args.boundary_order,args.basis_order,'dfull_basic_polynomials','DShape','dphi',[(1,0),(0,1)],rows)
    out += function(args.boundary_order,args.basis_order,'d2full_basic_polynomials','DShape','d2phi',[(2,0),(1,1),(0,2)],rows)
    print('\n'.join(out))

if __name__=='__main__': main()
