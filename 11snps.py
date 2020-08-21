from Bio import SeqIO
import re

#Estrutura de Dados Trie - Busca a maior sequência conservada de aminoácios em cada site já identificado
class trieNo(object):

    def __init__(self, chtr):
        self.chr = chtr
        self.seq = []
        self.fimpalavra = False
        self.contaSeq = 1
        self.contaAmino = 1
        self.porcenCon = 60

    #lê a sequencia de aminos no arquivo Fasta
    def leSequencias(self):
        listaSequencias = []
        for i in SeqIO.parse('seq-prot-ip3-R3.fas', 'fasta'):
            listaSequencias.append(list(i.seq))
        return listaSequencias


    #formata e adiciona os nomes das espécies à uma lista
    def leEspecies(self):
        listaEspecies = []
        for i in SeqIO.parse('seq-prot-ip3-R3.fas', 'fasta'):
            loc = re.search(r'type_3', i.id).end()
            sp = i.id[loc + 1:len(i.id)]
            sp = sp.split('_')
            sep = ' '
            item = sep.join(sp)
            listaEspecies.append(item)
        return listaEspecies


    def contaEspecies(self):
        listEsp = len(self.leEspecies())
        return listEsp



    # porcentagem a se considerar como site conservado - 60%
    def porcConservacao(self, numeroEspecies):
        return int(numeroEspecies * self.porcenCon/100)


    # elimina sequências com baixa conservação
    def locaisConserv(self, porConserv):
        listaSequencias = []
        for i in SeqIO.parse('seq-prot-ip3-R3.fas', 'fasta'):
            listaSequencias.append(i.seq)
        locaisConservados = []
        sitesConservados = []

        for j in range(len(listaSequencias[0])):
            linhas = []
            for k in range(len(listaSequencias)):
                if (listaSequencias[k][j] != '-'):
                    linhas.append(listaSequencias[k][j])

            if len(linhas) >= porConserv:
                sitesConservados.append(linhas.copy())
                locaisConservados.append(j)
            linhas.clear()
        return locaisConservados


    # agrupa os sites conservados consecutivos
    def agrupaSitesCons(self, listCon):
        listaTemp = []
        listaFim = []
        m = 0
        incr = 0
        tam = len(listCon)

        while incr < tam:
            chave = incr
            while listCon[incr] - listCon[chave] == m:
                listaTemp.append(listCon[incr])
                m = m + 1
                incr = incr + 1
                if incr == tam:
                    break
            listaFim.append(listaTemp.copy())
            listaTemp.clear()
            m = 0
        return listaFim


    # retorna uma lista com os sites conservados por espécie, uma linha por espécie
    def listAminosConservados(self, lisSitesCons, lisSeqAminos):
        listaTemp = []
        listaFim = []

        for i in range(len(lisSeqAminos)):
            for j in range(len(lisSeqAminos[i])):
                if j in lisSitesCons and lisSeqAminos[i][j] != '-':
                    listaTemp.append(lisSeqAminos[i][j])
            listaFim.append(listaTemp.copy())
            listaTemp.clear()
        return listaFim

    def maiorSequenciaConservada(self, lisCons, listSeq):
        maior = 0
        indice = 0
        listaFim = []
        for i in range(len(lisCons)):
            if len(lisCons[i]) > maior:
                maior = len(lisCons[i])
                indice = i

        inicio = lisCons[indice][0]
        final = lisCons[indice][maior-1]

        indNaoVazio = 0

        for j in range(len(listSeq)):
            if listSeq[j][inicio] != '-':
                indNaoVazio = j
                break
        listaFim.append(listSeq[indNaoVazio][inicio:final+1])
        sep = ''
        saida = sep.join(listaFim[0])


        return [maior,inicio, final, saida]

    def adicionaNo(self, palavra):
        node = self
        for caract in palavra:
            achou = False
            for item in node.seq:
                if item.chr == caract:
                    item.contaSeq = item.contaSeq+1
                    node = item
                    achou = True
                    break
            if not achou:
                novoNo = trieNo(caract)
                node.seq.append(novoNo)
                node = novoNo
        node.fimpalavra = True
        return [node.contaSeq, palavra]




    #pega os aminoacidos conservados por linha dentro da taxa de conservação estabelecida:
    def geraConserv(self, listaCon, listaSeq):
        listTemp = []
        listaAdd = []


        for item in listaCon:
            min = item[0]
            tam = len(item) -1
            max = item[tam]
            for p in range(len(listaSeq)):
                listTemp.append(listaSeq[p][min:max+1])

            listaAdd.append(listTemp.copy())
            listTemp.clear()


        #converte cada linha em uma palavra pra ser usada na Trie
        lTemp = []
        listaStr = []
        separator = ''
        for item in listaAdd:
            for j in item:
                strAm = separator.join(j)
                lTemp.append(strAm)
            listaStr.append(lTemp.copy())
            lTemp.clear()

        return listaStr



#gera as sequências mais encontradas em cada grupo de sites
    def geraPred(self, listStr):

        seqPred = ''
        listaSaidaPred = []

        for item in listStr:
            maior = 1
            for j in item:
                value = trieNo.adicionaNo(self, j)[0]
                nulo = j.find('-')
                if value > maior and nulo <= -1:
                    maior = value
                seqPred = trieNo.adicionaNo(self, j)[1]
            listaSaidaPred.append(seqPred)

        return listaSaidaPred

    def compara(raiz, listaAminos, seqPred, listEsp):

        result = []
        listPred = []
        listaConservada = []
        listEspAlter = []

        node = raiz
        raiz.adicionaNo(seqPred)


        if not raiz.seq:
            return False, 0

        esp = 0



        for value in listaAminos:
            sep = ''
            alterado = False
            x = 0
            while node.seq and x<len(value):
                node = node.seq[0]
                if node.chr == value[x]:
                    result.append(value[x])
                else:
                    result.append('[{}]'.format(value[x]))
                    alterado = True
                x = x+ 1
            if alterado:
                listPred.append(sep.join(result))
                listEspAlter.append(listEsp[esp])
            else:
                listaConservada.append(listEsp[esp])

            result.clear()
            node = raiz
            esp = esp+1
        return [listEspAlter, listPred, listaConservada]



#extrai as sequências diferentes das predominante
    def extraiSeqs(self, listEspecies, listaAminCons, listaPred):

        listaEspAlt = []
        listaAminosAlterados = []
        listaConservada = []
        listaEspConsTemp = []
        listaEspCons = []



        for i in range(len(listaPred)):
            for j in range(len(listaAminCons)):
                for k in range(len(listaAminCons[j])):
                    if listaPred[i] != listaAminCons[j][k]:
                        listaEspAlt.append(listEspecies[k])
                        listaAminosAlterados.append(listaAminCons[j][k])
                    else:
                        listaConservada.append(listaAminCons[j][k])
                        listaEspConsTemp.append(listEspecies[k])

                listaEspCons.append(listaEspConsTemp.copy())




        return [listaEspAlt, listaAminosAlterados, listaEspCons, listaConservada]


#identifica, nos sites alterados, qual amino foi mudado
    def altAminos(self, listSeq, listaCon, listaPred):

        for item in listaCon:
            min = item[0]
            tam = len(item) - 1
            max = item[tam]
            for p in range(len(listaSeq)):
                listTemp.append(listaSeq[p][min:max + 1])
            listaAdd.append(listTemp.copy())
            listTemp.clear()



        #for i in range(listaPred):
            #for j in range(listaPred[i]):







root = trieNo('+')

especiesList = root.leEspecies()
qtdeEspecies = root.contaEspecies()
porceCon = root.porcConservacao(qtdeEspecies)
locaisConservados = root.locaisConserv(porceCon)
seqAminos = root.leSequencias()
indices = root.agrupaSitesCons(locaisConservados)
conservados = root.geraConserv(indices, seqAminos)
predominantes = root.geraPred(conservados)



#cp = root.compara(conservados[0], predominantes[0], especiesList)




#extrai = root.extraiSeqs(especiesList, conservados, predominantes)




#alterados = root.altAminos(extrai[0], extrai[1], predominantes)


# gera arquivo de saida




quant = root.maiorSequenciaConservada(indices,seqAminos)[0]
inicio = root.maiorSequenciaConservada(indices, seqAminos)[1]
fim = root.maiorSequenciaConservada(indices, seqAminos)[2]
seqMaior = root.maiorSequenciaConservada(indices, seqAminos)[3]







with open("saida.txt", "w") as arq:
    text01 = '\nMaior sequência conservada\n'
    text02 = '\n - Quantidade de aminoácidos: {} \n'.format(quant)
    text03 = '\n - Sites da maior sequência conservada - site {} ao {}\n'.format(inicio + 1, fim + 1)
    text07 = '\nQuantidade de espécies existentes no alinhamento: {}\n'.format(qtdeEspecies)
    text04 = '\n - Sequência de aminoácidos da maior sequência conservada:\n'
    text05 = '\n - Porcentagem considerada para sequências conservadas: {}\n'.format(root.porcenCon)
    text06 = '\n - Quantidade mínima de espécies a considerar como região conservada\n'
    text08 = '\n(porcentagem X quant. espécies): {}\n'.format(porceCon)
    arq.write(text07)
    arq.write(text01)
    arq.write(text02)
    arq.write(text03)
    arq.write(text04)
    arq.write(seqMaior)
    arq.write('\n')
    arq.write(text05)
    arq.write(text06)
    arq.write(text08)
    arq.write('\n')
    arq.write('\n')



    arq.close()

# lista sites agrupados


# sequência predominante
# espécies que possuem a sequência predominante
# espeéciea alteradas com sua sequência
# amino alterado












"""
Fluxo:

- Lê sequências 
- Lê espécies
- Conta quantidade de espécies
- Define porcentagem de conservação a considerar - porConserv
- Calcula a maior sequência conservada - maiorSequenciaConservada(lisCons)
- Elimine da lista os sites com baixa conservação - locaisConserv(porConserv)
- Sequencia os sites conservados - agrupaSitesCons(listCon):
- Gera lista com os aminos para cada grupo de sites conservados sequnciados anteriormente e Concatena os aminos, 
  convertendo em palavras, para serem usados na Trie - - geraConserv(listaCon, listaSeq):
  
- Insere a lista com os nomes na função adicionaNo para busca da maior sequência - geraPred(listaStr)
- Compara lista gerada de aminos nos sites conservados com a sequência predominate - compara(listaAminos, seqPred)
- Gerar arquivo de saida


"""

















