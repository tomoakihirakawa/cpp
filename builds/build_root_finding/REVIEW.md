<div style="overflow: hidden;">
<p align="center">
  <img src="sample_aquarium.gif" style="object-fit: cover; width: 40%; height: 40%; margin-bottom: -20%;">
</p>
</div>

### はじめに

筋電図を使って魚の筋力分布を測定したところ，力は主に前部と中部の筋肉で発生しており，
多くの魚が持つ後部の細い尾柄はその力を主に後方へと伝達する役割を持っていることがわかっている．
魚の泳ぎは複雑で，様々なアプローチから研究されてきたが，
多くの場合，Lighthillの*細長い体の理論(Elongated Body Theory)*，またはその発展版が使われている．
([Lighthill 1969](https://www.annualreviews.org/doi/10.1146/annurev.fl.01.010169.002213)，[Lighthill 1971](https://royalsocietypublishing.org/doi/10.1098/rspb.1971.0085)，[Porez et al. 2014](https://journals.sagepub.com/doi/abs/10.1177/0278364914525811))

REFERENCE: [Yong Zhong et al. 2018](https://ieeexplore.ieee.org/document/8329488)

Wu 1961の理論は，無限高の板の波動運動に基づき，
一方，Lighthill 1960の理論は，線長い体の理論(Batchelor 1967 section 6.9を参照)の拡張としてelongated-body theory(EBT)に基づいている．

EBTの発展版として，Newman & Wu 1973やNewman1973がフィンを考慮した計算を行っている．
また，Lighthill 1970では，EBTをlarge-amplitude elongated-body theory(LAEBT)に拡張している．
LAEBTは，生体流体力学の分野で，魚の泳ぎの基礎となる理論として広く受け入れられている(Weihs 1972, Lighthill 1975).
ロボットに関しては，Boyer et al. 2008を参照．

しかし，数値解析結果との比較を通して，LAEBTは，実際の流れと大きく異なることがわかった(Wolfgang 1999, Triantafyllou et al. 2000)．

そこで，[Candelier et al. 2018](https://www.cambridge.org/core/product/identifier/S002211201000649X/type/journal_article)は，LAEBTに対応する実際の流れを調べ，さらにその３次元版を提案している．


REFERENCE: [Candelier et al. 2018](https://www.cambridge.org/core/product/identifier/S002211201000649X/type/journal_article)
