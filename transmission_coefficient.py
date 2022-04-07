# calculate transmission coefficient
import cmath
from units import *

def probe_ke(E, V): # 单位是hatree
    return cmath.sqrt(2*params['probe']['EM']*(params['probe']['fermi'] + E + V + params['probe']['workfunction']-params['superconductor']['Dirac_workfunction'])/Joule/Kilogram)/hbar_SI*a0_SI
def probe_kh(E, V):
    return -probe_ke(-E,V)
def vacuum_ke(E, V):
    return cmath.sqrt(2*params['probe']['EM']*( E + V/2 - params['superconductor']['Dirac_workfunction'])/Joule/Kilogram)/hbar_SI*a0_SI
def vacuum_kh(E, V):
    return -vacuum_ke(-E,V)
# omega的值为9.99E-6Ha
omega = (params['superconductor']['EM']/Kilogram)*((params['superconductor']['a']/Meter)**2)*((params['superconductor']['gap']/Joule)**2)/(hbar_SI**2)*Joule
def Sigma(E):
    return E**2 - 2*omega*params['superconductor']['Dirac_fermi'] + omega**2
def q1(E):
    return cmath.sqrt(2*params['superconductor']['EM']*(params['superconductor']['Dirac_fermi'] - omega + cmath.sqrt(Sigma(E)))/Joule/Kilogram)/hbar_SI*a0_SI
def q2(E):
    return cmath.sqrt(2*params['superconductor']['EM']*(params['superconductor']['Dirac_fermi'] - omega - cmath.sqrt(Sigma(E)))/Joule/Kilogram)/hbar_SI*a0_SI
# 波函数参数
def waveFunctionParams(E, type = 'electron', dir = 'up'):
    if (type == 'electron'):
        epsilon = -omega + cmath.sqrt(Sigma(E))
    elif (type == 'hole') :
        epsilon = -omega - cmath.sqrt(Sigma(E))
    else:
        raise Exception('type 只能是electron或者是hole'.format(type))
    if (dir == 'up'):
        if (Sigma(E) >= 0 and abs(E) >= abs(epsilon)): # 等价于--Sigma大于等于0, 可以判断出E的范围
            return cmath.sqrt((E+epsilon)/(2*abs(E)))
        elif (Sigma(E) >= 0 and abs(E) < abs(epsilon)): 
            return cmath.sqrt((E+epsilon)/(2*epsilon))
        else:
            return cmath.sqrt((E+epsilon)/(cmath.sqrt(2*omega*(params['superconductor']['Dirac_fermi']-E))+cmath.sqrt(2*omega*(params['superconductor']['Dirac_fermi']+E))))
    if (dir == 'down'):
        if (Sigma(E) >= 0 and abs(E) >= abs(epsilon)): 
            return cmath.sqrt((E-epsilon)/(2*abs(E)))
        elif (Sigma(E) >= 0 and abs(E) < abs(epsilon)): 
            return cmath.sqrt((E-epsilon)/(2*epsilon))
        else:
            return cmath.sqrt((E-epsilon)/(cmath.sqrt(2*omega*(params['superconductor']['Dirac_fermi']-E))+cmath.sqrt(2*omega*(params['superconductor']['Dirac_fermi']+E))))
    else:
        raise Exception('dir 只能是up或者是down'.format(dir))
# 总结：定义的函数probe_ke、probe_kh、vacuum_ke、vacuum_kh是E和V的函数
#       Sigma、q1、q2、waveFunctionParams是E的函数
def transmissionParams(E,Ee,Eh, V, d, type='ref', state=1): # 默认值是态1的反射系数
    # 分母部分
    deno = waveFunctionParams(E,'electron','up')*waveFunctionParams(E,'hole','down')*cmath.exp(1.0j*probe_ke(Ee,V)*d/2)*(
        (probe_ke(Ee,V)+vacuum_ke(Ee,V))*(vacuum_ke(Ee,V)+q1(E))*cmath.exp(-1.0j*vacuum_ke(Ee,V)*d) + (probe_ke(Ee,V)-vacuum_ke(Ee,V))*(vacuum_ke(Ee,V)-q1(E))*cmath.exp(1j*vacuum_ke(Ee,V)*d)
    )*(
        (probe_kh(Eh,V)+vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)+q2(E))*cmath.exp(-1.0j*vacuum_kh(Eh,V)*d) + (probe_kh(Eh,V)-vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)-q2(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
    )-\
    waveFunctionParams(E,'electron','down')*waveFunctionParams(E,'hole','up')*cmath.exp(1j*probe_ke(Ee,V)*d/2)*(
        (probe_ke(Ee,V)+vacuum_ke(Ee,V))*(vacuum_ke(Ee,V)+q2(E))*cmath.exp(-1j*vacuum_ke(Ee,V)*d) + (probe_ke(Ee,V)-vacuum_ke(Ee,V))*(vacuum_ke(Ee,V)-q2(E))*cmath.exp(1j*vacuum_ke(Ee,V)*d)
    )*(
        (probe_kh(Eh,V)+vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)+q1(E))*cmath.exp(-1j*vacuum_kh(Eh,V)*d) + (probe_kh(Eh,V)-vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)-q1(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
    )
    #透射系数和反射系数
    if (type == "trans"):
        deno = deno/cmath.exp(1j*(q1(E)+q2(E))*d/2)
        if (state == 1):
            numerator = 4*probe_ke(Ee,V)*vacuum_ke(Ee,V)*waveFunctionParams(E,'hole','down')*cmath.exp(1j*q2(E)*d/2)*(
                (probe_kh(Eh,V)+vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)+q2(E))*cmath.exp(-1j*vacuum_kh(Eh,V)*d) + (probe_kh(Eh,V)-vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)-q2(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
            )
            return numerator/deno
        elif (state == 2):
            numerator = 4*probe_ke(Ee,V)*vacuum_ke(Ee,V)*(-waveFunctionParams(E,'electron','down'))*cmath.exp(1j*q1(E)*d/2)*(
                (probe_kh(Eh,V)+vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)+q1(E))*cmath.exp(-1j*vacuum_kh(Eh,V)*d) + (probe_kh(Eh,V)-vacuum_kh(Eh,V))*(vacuum_kh(Eh,V)-q1(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
            )
            return numerator/deno
        else:
            raise Exception('state 只能为1或2'.format(state))
    elif (type == "ref"):
        if (state == 1):
            numerator = waveFunctionParams(E,'electron','up')*waveFunctionParams(E,'hole','down')*cmath.exp(-1j*probe_ke(Ee,V)*d/2)*(
                (
                    (probe_ke(Ee,V)-vacuum_ke(Ee,V)*(vacuum_ke(Ee,V)+q1(E)))*cmath.exp(-1j*vacuum_ke(Ee,V)*d)+(probe_ke(Ee,V)+vacuum_ke(Ee,V)*(vacuum_ke(Ee,V)-q1(E)))*cmath.exp(1j*vacuum_ke(Ee,V)*d)
                )*(
                    (probe_kh(Eh,V)+vacuum_kh(Eh,V)*(vacuum_kh(Eh,V)+q2(E)))*cmath.exp(-1j*vacuum_kh(Eh,V)*d)+(probe_kh(Eh,V)-vacuum_kh(Eh,V)*(vacuum_kh(Eh,V)-q2(E)))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
                )
            )-waveFunctionParams(E,'electron','down')*waveFunctionParams(E,'hole','up')*cmath.exp(-1j*probe_ke(Ee,V)*d/2)*(
                (
                    (probe_ke(Ee,V)-vacuum_ke(Ee,V)*(vacuum_ke(Ee,V)+q2(E)))*cmath.exp(-1j*vacuum_ke(Ee,V)*d)+(probe_ke(Ee,V)+vacuum_ke(Ee,V)*(vacuum_ke(Ee,V)-q2(E)))*cmath.exp(1j*vacuum_ke(Ee,V)*d)
                )*(
                    (probe_kh(Eh,V)+vacuum_kh(Eh,V)*(vacuum_kh(Eh,V)+q1(E)))*cmath.exp(-1j*vacuum_kh(Eh,V)*d)+(probe_kh(Eh,V)-vacuum_kh(Eh,V)*(vacuum_kh(Eh,V)-q1(E)))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
                )
            )
            return numerator/deno
        elif(state == 2):
            numerator = (2*vacuum_ke(Ee,V)/vacuum_kh(Ee,V))*(
                probe_ke(Ee,V)*waveFunctionParams(E,'electron','down')*waveFunctionParams(E,'hole','down')*cmath.exp(-1j*probe_kh(Eh,V)*d/2)
            )*(
                (
                    (probe_kh(Eh,V)+q1(E))*cmath.exp(-1j*vacuum_kh(Eh,V)*d)+(vacuum_kh(Eh,V)-q1(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
                )*(
                    (probe_kh(Eh,V)+vacuum_kh(Eh,V))*(probe_kh(Eh,V)+q2(E))*cmath.exp(-1j*vacuum_kh(Eh,V)*d)+(probe_kh(Eh,V)-vacuum_kh(Eh,V))*(probe_kh(Eh,V)-q2(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
                )
            )-(
                probe_ke(Ee,V)*waveFunctionParams(E,'electron','down')*waveFunctionParams(E,'hole','down')*cmath.exp(-1j*probe_kh(Eh,V)*d/2)
            )*(
                (
                    (probe_kh(Eh,V)+q2(E))*cmath.exp(-1j*vacuum_kh(Eh,V)*d)+(vacuum_kh(Eh,V)-q2(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
                )*(
                    (probe_kh(Eh,V)+vacuum_kh(Eh,V))*(probe_kh(Eh,V)+q1(E))*cmath.exp(-1j*vacuum_kh(Eh,V)*d)+(probe_kh(Eh,V)-vacuum_kh(Eh,V))*(probe_kh(Eh,V)-q1(E))*cmath.exp(1j*vacuum_kh(Eh,V)*d)
                )
            )
            return numerator/deno
        else:
            raise Exception('state 只能为1或2'.format(state))
    else:
        raise Exception('type 只能是反射ref或者是透射trans'.format(type))
# 测试：
E = 0*eV
Ee = 0*eV
Eh = 0*eV
V = 0*eV
print(V*27)
d = 0.5e-10*Meter
print('ke=', probe_ke(Ee,V)*Angstrom)
print('kh=', probe_kh(Eh,V)*Angstrom)
print('ve=', vacuum_ke(Ee,V)*Angstrom)
print('vh=', vacuum_kh(Eh,V)*Angstrom)
print('q1=',q1(E)*Angstrom)
print('q2=',q2(E)*Angstrom)
print('u01=', waveFunctionParams(E,'electron','up'))
print('u02=', waveFunctionParams(E,'hole','up'))
print('v01=', waveFunctionParams(E,'electron','down'))
print('v02=', waveFunctionParams(E,'hole','down'))
print('re1=', transmissionParams(E,Ee,Eh,V,d,'ref', 1),'and |re1|=',abs(transmissionParams(E,Ee,Eh,V,d,'ref', 1)))
print('re2=', transmissionParams(E,Ee,Eh,V,d,'ref', 2),'and |re2|=',abs(transmissionParams(E,Ee,Eh,V,d,'ref', 2)))
print('te1=', transmissionParams(E,Ee,Eh,V,d,'trans', 1),'and |te1|=',abs(transmissionParams(E,Ee,Eh,V,d,'trans', 1)))
print('te2=', transmissionParams(E,Ee,Eh,V,d,'trans', 2),'and |te2|=',abs(transmissionParams(E,Ee,Eh,V,d,'trans', 2)))