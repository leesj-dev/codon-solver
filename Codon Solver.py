# 필요한 패키지 불러오기
from itertools import product
from operator import itemgetter
import re
import numpy

# 아미노산(한글) to 아미노산(one-letter code)
korean_dict = {"페닐알라닌": "f", "류신": "l", "아이소류신": "i", "메싸이오닌": "m", "발린": "v", "세린": "s", "프롤린": "p", "트레오닌": "t", "알라닌": "a", "타이로신": "y", "종결코돈": "*", "히스티딘": "h", "글루타민": "q", "아스파라진": "n", "라이신": "k", "아스파트산": "d", "글루탐산": "e", "시스테인": "c", "트립토판": "w", "아르지닌": "r", "글라이신": "g"}

# 아미노산(one-letter code) to 아미노산(한글): key와 value를 서로 바꿈
oneletter_dict = dict((value, key) for key, value in korean_dict.items())

# 아미노산 to 코돈 변환 dictionary
amino_dict = {
    "m": ["AUG"],
    "w": ["UGG"],
    "f": ["UUU", "UUC"],
    "h": ["CAU", "CAC"],
    "q": ["CAA", "CAG"],
    "n": ["AAU", "AAC"],
    "k": ["AAA", "AAG"],
    "d": ["GAU", "GAC"],
    "e": ["GAA", "GAG"],
    "c": ["UGU", "UGC"],
    "y": ["UAU", "UAC"],
    "i": ["AUU", "AUC", "AUA"],
    "*": ["UAA", "UAG", "UGA"],
    "v": ["GUU", "GUC", "GUA", "GUG"],
    "p": ["CCU", "CCC", "CCA", "CCG"],
    "t": ["ACU", "ACC", "ACA", "ACG"],
    "a": ["GCU", "GCC", "GCA", "GCG"],
    "g": ["GGU", "GGC", "GGA", "GGG"],
    "s": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "l": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "r": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"] 
}

# 코돈 to 아미노산 변환 dictionary
triplet_dict = {
    "AUG" : "m",
    "UGG" : "w",
    "UUU" : "f", "UUC" : "f",
    "CAU" : "h", "CAC" : "h",
    "CAA" : "q", "CAG" : "q",
    "AAU" : "n", "AAC" : "n",
    "AAA" : "k", "AAG" : "k",
    "GAU" : "d", "GAC" : "d",
    "GAA" : "e", "GAG" : "e",
    "UGU" : "c", "UGC" : "c",
    "UAU" : "y", "UAC" : "y",
    "AUU" : "i", "AUC" : "i", "AUA" : "i",
    "UAA" : "*", "UAG" : "*", "UGA" : "*",
    "GUU" : "v", "GUC" : "v", "GUA" : "v", "GUG" : "v",
    "CCU" : "p", "CCC" : "p", "CCA" : "p", "CCG" : "p",
    "ACU" : "t", "ACC" : "t", "ACA" : "t", "ACG" : "t",
    "GCU" : "a", "GCC" : "a", "GCA" : "a", "GCG" : "a",
    "GGU" : "g", "GGC" : "g", "GGA" : "g", "GGG" : "g",
    "UCU" : "s", "UCC" : "s", "UCA" : "s", "UCG" : "s", "AGU" : "s", "AGC" : "s",
    "UUA" : "l", "UUG" : "l", "CUU" : "l", "CUC" : "l", "CUA" : "l", "CUG" : "l",
    "CGU" : "r", "CGC" : "r", "CGA" : "r", "CGG" : "r", "AGA" : "r", "AGG" : "r"
    }

# 잉여 문자 제거 dictionary
leftover_dict1 = {"'": "", ",": "", "[": "", "]": "", " ": ""}
leftover_dict2 = {"'": "", '"': "", "[": "", "]": "", "(, ": "(", ", )": ")"}

# string 안에 들어 있는 여러 개의 substring을 replace하는 함수
# rep이라는 dictionary 안에 있는 각각의 key를 text에서 찾은 후 그 key에 해당하는 value로 대체함
def replace_all(text, rep):
    rep = dict((re.escape(k), v) for k, v in rep.items())
    pattern = re.compile("|".join(rep.keys()))
    text = pattern.sub(lambda m: rep[re.escape(m.group(0))], text)
    return text

# RNA 서열을 3개 단위로 쪼개어 list를 만드는 함수
def Split(RNA):
    split = []
    for i in range(0, len(RNA), 3):
        split.append(RNA[i : i + 3])
    return split

# RNA 서열을 아미노산으로 번역해주는 함수
def convert(RNA):
    split = Split(RNA)                              # RNA 쪼개기
    triplet = str(split)                            # string형으로 변환
    amino = replace_all(triplet, triplet_dict)      # 아미노산으로 변환
    amino = replace_all(amino, leftover_dict1)      # 잉여 문자 제거
    return amino

# DNA로 역전사
def tran(rna):
    dna = rna.replace("C", "c")
    dna = dna.replace("G", "C")
    dna = dna.replace("c", "G")
    dna = dna.replace("A", "T")
    dna = dna.replace("U", "A")
    return dna

# 폴리펩타이드(one-letter code)를 폴리펩타이드(한글)로 바꿔 줌
def protein_process(protein):
    protein = '-'.join(protein)
    protein = replace_all(protein, oneletter_dict)
    return protein

# X로부터 합성된 폴리펩타이드 서열 (편의를 위해 아미노산을 one-letter code로 적음)
xProtein_ko = "메싸이오닌-발린-라이신-히스티딘-아스파트산-류신-세린-아르지닌"

# xProtein(한글)을 xProtein(one-letter code)로 바꾸고 잉여 문자 제거
xProtein = replace_all(xProtein_ko, korean_dict)
xProtein = xProtein.replace("-", "")

# amino_dict에서 xProtein에 해당하는 코돈을 찿아 경우의 수를 xRNA_list 안에 모두 나열
xRNA_lst = list(map("".join, product(*itemgetter(*xProtein)(amino_dict))))

# 5개의 뉴클레오타이드, 14개의 수소 결합일 때 A-T 결합 수(a)와 G-C 결합 수(b) 연립방정식 풀기 (행렬 이용)
A = numpy.array([[1, 1], [2, 3]])
B = numpy.array([5, 14])
C = numpy.linalg.solve(A, B)
a = int(C[0])
b = int(C[1])

# 앞으로 채워넣을 빈 list 자료형 생성
DNA_list = []           # 문제의 조건을 만족하는 xDNA, yDNA, zDNA 순서쌍을 담을 list
RNA_list = []           # 문제의 조건을 만족하는 xRNA, yRNA, zRNA 순서쌍을 담을 list
Protein_list = []       # 문제의 조건을 만족하는 X, Y, Z 순서쌍을 담을 list
dcodon_list = []        # 문제의 조건을 만족하는 X의 '(가)아스파트산'에 해당하는 코돈들
num3_list = []          # ㄷ 선지 해결을 위한 list

for xRNA in xRNA_lst:                                                                                               # xRNA_lst 안에 들어있는 xRNA 후보군 각각에 대하여
    for j in range(0, len(xRNA) - 5):                                                                               # j = 첫번째 글자부터 (xRNA의 길이 - 4)번째 글자까지 반복           
        part = xRNA[j : j + 5]                                                                                      # xRNA의 j번째부터 j+5번째 글짜까지 5개 글자만 추출
        bond2 = part.count("A") + part.count("U")                                                                   # A-T 결합 수 세기
        bond3 = part.count("G") + part.count("C")                                                                   # G-C 결합 수 세기
        if bond2 == a and bond3 == b:                                                                               # 잘라낸 부분이 a개의 A-T 결합 수와 b개의 G-C 결합 수를 포함한다면
            yRNA = xRNA.replace(xRNA[j : j + 5], "")                                                                # yRNA = xRNA에서 연속된 그 5개의 뉴클레오타이드 결실
            yProtein = convert(yRNA)                                                                                # Y 추출
            yProtein = yProtein.split("*", 1)[0]                                                                    # Y의 종결코돈에 해당하는 "*" 및 그 이후 아미노산 제거
            if (yProtein[0] == "m" and len(yProtein) == 4):                                                         # Y가 메싸이오닌에서 시작하고 아미노산 4개로 이루어져 있다면
                cnt = 0                                                                                             # cnt 변수 초기화
                if len(set(yProtein)) == len(yProtein):                                                             # Y에 있는 아미노산이 모두 서로 다르다면
                    for ltr in yProtein:                                                                            # Y 중 X에도 있는 아미노산 개수(cnt) 세기
                        if ltr in xProtein:
                            cnt = cnt + 1
                    if cnt == 3:                                                                                    # 만약 Y 중 3개의 아미노산만 X에 있다면
                        for k in range(0, len(xRNA)):                                                               # xRNA의 첫번째 글자부터 마지막 글자까지 반복
                            zRNA = (xRNA[:k] + "UU" + xRNA[k:])                                                     # xRNA에 UU를 삽입하여 zRNA를 만들기 (xDNA에 AA가 삽입되었으므로)
                            zProtein = convert(zRNA)                                                                # Z 추출
                            zProtein = zProtein.split("*", 1)[0]                                                    # Z의 종결코돈에 해당하는 "*" 및 그 이후 아미노산 제거
                            if (zProtein[0] == "m" and len(zProtein) < 4):                                          # Z가 메싸이오닌에서 시작하고 Z가 Y보다 적은 수의 아미노산으로 이루어져 있다면
                                RNA_list.append(str(["(" + xRNA, yRNA, zRNA + ")"]))                                # RNA_list에 해당 xRNA, yRNA, zRNA 순서쌍을 추가

                                xProtein_ko = protein_process(xProtein)                                             # 폴리펩타이드를 한글로 바꿈
                                yProtein_ko = protein_process(yProtein)
                                zProtein_ko = protein_process(zProtein)
                                Protein_list.append(str(["(" + xProtein_ko, yProtein_ko, zProtein_ko + ")"]))       # Protein_list에 한글로 바꾼 X, Y, Z 순서쌍을 추가

                                xDNA = tran(xRNA)                                                                   # xDNA, yDNA, zDNA 각각을 RNA로부터 역전사
                                yDNA = tran(yRNA)
                                zDNA = tran(zRNA)
                                DNA_list.append(str(["(" + xDNA, yDNA, zDNA + ")"]))                                # DNA_list에 해당 xDNA, yDNA, zDNA 순서쌍을 추가
                                print(xDNA, xProtein, yDNA, yProtein, zDNA, zProtein)                               # x~z의 DNA와 폴리펩타이드 경우의 수를 출력 (가독성을 위해 one-letter code 아미노산 사용)

                                # ㄱ 해결: 문제 조건을 만족하는 모든 경우의 수에 대하여 (가)의 코돈이 GAC만 있는지 확인하기 위해 list에 후보군을 계속 추가하는 방식을 사용
                                xSplit = Split(xRNA)                                                                # xRNA를 3개 단위로 쪼개기
                                for item in xSplit:                                                                 # xRNA 각각의 코돈에 대하여
                                    if convert(item) == "d":                                                        # 그 코돈으로부터 번역된 아미노산이 아스파트산이면
                                        dcodon_list.append(item)                                                    # dcodon_list에 그 코돈 항목을 추가

                                # ㄷ 해결: 문제 조건을 만족하는 모든 경우의 수에 대하여 yRNA와 zRNA의 종결코돈이 동일한지에 관한 bool값을 list에 계속 추가하는 방식을 사용 
                                ySplit = Split(yRNA)                                                                # yRNA를 3개 단위로 쪼개기
                                for item in ySplit:                                                                 # yRNA 각각의 코돈에 대하여
                                    if convert(item) == "*":                                                        # 그 코돈이 종결코돈이라면
                                        ycodon = item                                                               # ycodon을 그 코돈으로 정의

                                zSplit = Split(zRNA)                                                                # zRNA를 3개 단위로 쪼개기
                                for item in zSplit:                                                                 # zRNA 각각의 코돈에 대하여
                                    if convert(item) == "*":                                                        # 그 코돈이 종결코돈이라면
                                        zcodon = item                                                               # zcodon을 그 코돈으로 정의

                                if ycodon == zcodon:                                                                # 만약 yRNA와 zRNA의 종결코돈이 동일하다면
                                    num3 = True                     
                                else:
                                    num3 = False                    

                                num3_list.append(num3)                                                              # num3_list에 num3 값을 추가
print("\n")

total_list = [DNA_list, RNA_list, Protein_list, dcodon_list, num3_list]     # 5개의 list로 구성된 새로운 list 생성
for i in range(len(total_list)):                                            # list 5개 각각에 대하여
    total_list[i] = str(list(dict.fromkeys(total_list[i])))                 # list의 중복값을 제거하기 위해, list에서 dictionary를 만들어 그 key를 다시 list로 바꿈
    total_list[i] = replace_all(total_list[i], leftover_dict2)              # 잉여 글자 제거
    print(total_list[i])                                                    # list 출력
    print("\n")                                                             # 한 칸 띄우기

# ㄱ
if total_list[3] == "GAC":      # dcodon_list가 GAC만 가지고 있다면
    No1 = True
else:
    No1 = False

# ㄴ
if len(zProtein) == 3:          # Z가 아미노산 3개로 이루어져 있다면
    No2 = True
else:
    No2 = False

# ㄷ
if total_list[4] == "True":     # num3_list가 참값만 가지고 있다면
    No3 = True
else:
    No3 = False

# 최종 답 출력
print("ㄱ.", No1, ", ", "ㄴ.", No2, ", ", "ㄷ.", No3)