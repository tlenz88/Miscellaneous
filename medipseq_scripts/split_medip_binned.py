import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep="\t", header=None)
b1s = []
b2s = []
b3s = []
b4s = []
b5s = []
b6s = []
b7s = []
b8s = []
b9s = []
b10s = []
b11s = []
b12s = []
b13s = []
b14s = []
b15s = []
b1e = []
b2e = []
b3e = []
b4e = []
b5e = []
b6e = []
b7e = []
b8e = []
b9e = []
b10e = []
b11e = []
b12e = []
b13e = []
b14e = []
b15e = []
for i in df.itertuples():
	bs = (i[5] - i[4]) / 5
	if i[7] == "+":
		b1s.append(i[4]-500)
		b2s.append(i[4]-400)
		b3s.append(i[4]-300)
		b4s.append(i[4]-200)
		b5s.append(i[4]-100)
		b6s.append(i[4])
		b7s.append(int(round(i[4]+bs)))
		b8s.append(int(round(i[4]+2*bs)))
		b9s.append(int(round(i[4]+3*bs)))
		b10s.append(int(round(i[4]+4*bs)))
		b11s.append(i[5]+1)
		b12s.append(i[5]+101)
		b13s.append(i[5]+201)
		b14s.append(i[5]+301)
		b15s.append(i[5]+401)

		b1e.append(i[4]-401)
		b2e.append(i[4]-301)
		b3e.append(i[4]-201)
		b4e.append(i[4]-101)
		b5e.append(i[4]-1)
		b6e.append(int(round(i[4]+bs-1)))
		b7e.append(int(round(i[4]+2*bs-1)))
		b8e.append(int(round(i[4]+3*bs-1)))
		b9e.append(int(round(i[4]+4*bs-1)))
		b10e.append(i[5])
		b11e.append(i[5]+100)
		b12e.append(i[5]+200)
		b13e.append(i[5]+300)
		b14e.append(i[5]+400)
		b15e.append(i[5]+500)

	if i[7] == "-":
		b1s.append(i[5]+401)
		b2s.append(i[5]+301)
		b3s.append(i[5]+201)
		b4s.append(i[5]+101)
		b5s.append(i[5]+1)
		b6s.append(int(round(i[5]-bs+1)))
		b7s.append(int(round(i[5]-2*bs+1)))
		b8s.append(int(round(i[5]-3*bs+1)))
		b9s.append(int(round(i[5]-4*bs+1)))
		b10s.append(i[4])
		b11s.append(i[4]-100)
		b12s.append(i[4]-200)
		b13s.append(i[4]-300)
		b14s.append(i[4]-400)
		b15s.append(i[4]-500)

		b1e.append(i[5]+500)
		b2e.append(i[5]+400)
		b3e.append(i[5]+300)
		b4e.append(i[5]+200)
		b5e.append(i[5]+100)
		b6e.append(i[5])
		b7e.append(int(round(i[5]-bs)))
		b8e.append(int(round(i[5]-2*bs)))
		b9e.append(int(round(i[5]-3*bs)))
		b10e.append(int(round(i[5]-4*bs)))
		b11e.append(i[4]-1)
		b12e.append(i[4]-101)
		b13e.append(i[4]-201)
		b14e.append(i[4]-301)
		b15e.append(i[4]-401)

starts = [b1s, b2s, b3s, b4s, b5s, b6s, b7s, b8s, b9s, b10s, b11s, b12s, b13s, b14s, b15s]
ends = [b1e, b2e, b3e, b4e, b5e, b6e, b7e, b8e, b9e, b10e, b11e, b12e, b13e, b14e, b15e]

df1 = df.copy()
df1[3] = b1s
df1[4] = b1e

df2 = df.copy()
df2[3] = b2s
df2[4] = b2e

df3 = df.copy()
df3[3] = b3s
df3[4] = b3e

df4 = df.copy()
df4[3] = b4s
df4[4] = b4e

df5 = df.copy()
df5[3] = b5s
df5[4] = b5e

df6 = df.copy()
df6[3] = b6s
df6[4] = b6e

df7 = df.copy()
df7[3] = b7s
df7[4] = b7e

df8 = df.copy()
df8[3] = b8s
df8[4] = b8e

df9 = df.copy()
df9[3] = b9s
df9[4] = b9e

df10 = df.copy()
df10[3] = b10s
df10[4] = b10e

df11 = df.copy()
df11[3] = b11s
df11[4] = b11e

df12 = df.copy()
df12[3] = b12s
df12[4] = b12e

df13 = df.copy()
df13[3] = b13s
df13[4] = b13e

df14 = df.copy()
df14[3] = b14s
df14[4] = b14e

df15 = df.copy()
df15[3] = b15s
df15[4] = b15e

df1.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin1.gff", sep="\t", index=None, header=None)
df2.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin2.gff", sep="\t", index=None, header=None)
df3.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin3.gff", sep="\t", index=None, header=None)
df4.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin4.gff", sep="\t", index=None, header=None)
df5.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin5.gff", sep="\t", index=None, header=None)
df6.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin6.gff", sep="\t", index=None, header=None)
df7.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin7.gff", sep="\t", index=None, header=None)
df8.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin8.gff", sep="\t", index=None, header=None)
df9.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin9.gff", sep="\t", index=None, header=None)
df10.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin10.gff", sep="\t", index=None, header=None)
df11.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin11.gff", sep="\t", index=None, header=None)
df12.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin12.gff", sep="\t", index=None, header=None)
df13.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin13.gff", sep="\t", index=None, header=None)
df14.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin14.gff", sep="\t", index=None, header=None)
df15.to_csv("/mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_genes_bin15.gff", sep="\t", index=None, header=None)