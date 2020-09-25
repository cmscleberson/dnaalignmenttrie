from Bio import SeqIO
import re




from dnaTrie.principal import trieNo

class compilado:


    root = trieNo('+')

    especiesList = root.leEspecies()
    qtdeEspecies = root.contaEspecies()
    porceCon = root.porcConservacao(qtdeEspecies)
    locaisConservados = root.locaisConserv(porceCon)
    seqAminos = root.leSequencias()
    indices = root.agrupaSitesCons(locaisConservados)
    conservados = root.geraConserv(indices, seqAminos)
    predominantes = root.geraPred(conservados)


    # gera arquivo de saida

    arquivo = "saidaComp - Tipo"+str(root.tipo)+".txt"
    with open(arquivo, "w") as arq:

        quant = root.maiorSequenciaConservada(indices, seqAminos)[0]
        inicio = root.maiorSequenciaConservada(indices, seqAminos)[1]
        fim = root.maiorSequenciaConservada(indices, seqAminos)[2]
        seqMaior = root.maiorSequenciaConservada(indices, seqAminos)[3]

        text = '\nAlinhamentos Para o Tipo:\n{}\n'.format(root.tipo)
        text01 = '\nMaior sequência conservada\n'
        text02 = '\n - Quantidade de aminoácidos: {} \n'.format(quant)
        text03 = '\n - Sites da maior sequência conservada - site {} ao {}\n'.format(inicio + 1, fim + 1)
        text07 = '\nQuantidade de espécies existentes no alinhamento: {}\n'.format(qtdeEspecies)
        text04 = '\n - Sequência de aminoácidos da maior sequência conservada:\n'
        text05 = '\n - Porcentagem considerada para sequências conservadas: {}\n'.format(root.porcenCon)
        text06 = '\n - Quantidade mínima de espécies a considerar como região conservada\n'
        text08 = '\n(porcentagem X quant. espécies): {}\n'.format(porceCon)

        arq.write(text)
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
        arq.write('\nOrdem das espécies mais conservadas:\n')

        # parâmtros - listPred, listAminos, listEsp:
        for i in range(len(indices)):
            root.compara(predominantes[i], conservados[i], especiesList)

        # return [vetDes, vetOrd, vetEspDes, vetEspOrd]
        esp = root.geraRankOrd(especiesList)[3]
        alt = root.geraRankOrd(especiesList)[1]

        for i in range(len(esp)):
            arq.write('Ord.: {} - Espécie: {} - Alterações nos sites conservados: {}'.format(i+1, esp[i],alt[i]))
            arq.write('\n')
        arq.write('\n')
        arq.write('--------------------------------------------------------------------------------------------------------------------------')
        arq.write('\n\n')

        arq.close()


    print(root.geraRankOrd(especiesList)[2])
    print(root.geraRankOrd(especiesList)[3])














