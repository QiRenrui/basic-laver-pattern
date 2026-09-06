# A 30-Year-Old Problem on Laver Tables

Laver tables give rise to computable functions that can be implemented by remarkably short programs. One such implementation is:

```python
def l(n, a, b):
    return b if a == 1 << n else a + 1 if b < 2 else l(n, l(n, a, b - 1), a + 1)


def F(n, m=1):
    k = 1
    for _ in range(n - 1):
        k = l(m, 1, k)
    return m - 1 if k - (1 << m) else F(n, m + 1)
```

Randall Dougherty gave a lower bound for the growth of the relevant Laver-table function in his paper [*Critical points in an algebra of elementary embeddings*](https://arxiv.org/abs/math/9205202). Although enormous by ordinary standards, his lower bound is still below many of the large numbers commonly discussed in googology. For more than thirty years, no substantial improvement to this lower-bound method was known.

The situation is striking because the finite combinatorial objects involved have an elementary description, while the only known proof of the totality of the associated function uses extraordinarily strong [rank-into-rank large-cardinal axioms](https://neugierde.github.io/cantors-attic/Rank_into_rank).

The proof-theoretic strength needed to prove the totality of a computable function constrains how fast that function can grow. For example, if the totality of a computable function can be proved in Peano Arithmetic (PA), then its growth is eventually bounded by an appropriate level of the fast-growing hierarchy below $\varepsilon_0$. Dougherty therefore asked in 1992 whether one could prove a lower bound strong enough to show that the Laver-table function eventually dominates the PA-provably total computable functions. Equivalently, this would show that its totality is not provable in PA.

This problem became well known among people interested in Laver tables and their connections with set theory. Dougherty and Thomas Jech worked extensively on the surrounding questions, and Patrick Dehornoy repeatedly discussed Laver tables and their unusual dependence on strong set-theoretic principles. The subject has also appeared in expository writing, including John Baez's [blog post](https://johncarlosbaez.wordpress.com/2016/05/06/shelves-and-the-infinite/) and a [MathOverflow discussion](https://mathoverflow.net/q/206291).

Dougherty participated in the IMO three times and won one gold and two silver medals. He later collaborated with Matthew Foreman on the Banach--Tarski paradox, resolving a problem about Baire-measurable pieces that had been open for decades. Jech is one of the central figures in modern set theory, and his book *Set Theory* is a standard reference. Dehornoy made fundamental contributions to braid groups and self-distributive algebra, where Laver tables play an important role. This gives some sense of the mathematical background surrounding the problem.

My paper resolves Dougherty's question: the function $F$ arising from the critical points in the Laver-table setting eventually dominates every computable function whose totality is provable in PA. This was one of two problems open for roughly thirty years that I worked on during my PhD.

The route to the result began, somewhat unexpectedly, with googology, and eventually returned to googology.

## From googology to Laver tables

I had been interested in googology---the study of extremely large finite numbers and the ordinal notations used to organize fast-growing functions---for years before I seriously studied Laver tables. The existence of these functions was already familiar in the googology community, but the known difficulty of the problem had discouraged me from trying to improve Dougherty's bounds.

At one point I noticed that the ordinal notation system [BMS (Bashicu Matrix System)](https://googology.miraheze.org/wiki/Bashicu_matrix_system) could be viewed, roughly speaking, as a version of [patterns of resemblance (PoR)](https://www.sciencedirect.com/science/article/pii/S0168007200000403) with downward reflection but without the corresponding upward reflection. I proved its well-foundedness using $\Sigma_n$-elementary embeddings.

A major theme of googology is the search for stronger recursive ordinal notation systems. This led me to explore large-cardinal principles that might naturally generate such notations through structures analogous to $\Sigma_n$-elementary embeddings. I tried several possibilities, including ineffable and Woodin cardinals, without finding a satisfactory construction. Rank-into-rank axioms were different: they naturally produce the finite combinatorial structures known as Laver tables.

The relevant Laver-table function can be understood through counting certain critical points. The critical points themselves are ordinals and hence are naturally well ordered. However, the order type of the critical points directly used in the definition is only $\omega$, so they do not immediately provide a strong ordinal notation system. I began looking for ways to enlarge the collection of critical-point configurations so that the resulting order type would become a much larger recursive ordinal.

This led to a recursive notation system that I called **Laver table yarn (LTY)**: it can be thought of as an interweaving of Laver-table structures.

I did some initial exploration of LTY and discussed it with friends in the Chinese googology community. During this stage, [HypCos](https://github.com/hypcos) helped me with computations. As I tried to estimate a lower bound for the strength of LTY, I realized that I first needed a much better understanding of lower bounds for the original Laver-table function. At that point I did not expect to solve Dougherty's problem; I only hoped that the known estimates would already give a useful lower bound for LTY.

## Basic Laver patterns

I began systematically studying lower bounds for Laver tables. My scratch paper, and even the drawing tool in my phone's notes app while I was riding buses, filled up with diagrams of circles representing critical-point configurations.

<!-- TODO: restore the original screenshot of hand-drawn Laver patterns. -->

Gradually I found a route that gave much stronger bounds than Dougherty's general estimate. At an intermediate stage I already had enough to reach the $\varepsilon_0$ range, which would settle the PA-independence question. From a googology perspective, however, I wanted a much more precise analysis.

I therefore worked through Dougherty's rather ad hoc constructions in detail. They are not easy to read, but after understanding them I was able to modify and extend them substantially.

The key combinatorial structure I introduced is the **basic Laver pattern (BLP)**. A BLP supports operations analogous to those used with recursive ordinals: one can distinguish zero, successor, and limit behavior, and define operations corresponding to fundamental sequences and predecessors. Certain lower bounds for the Laver-table function then arise from fast-growing hierarchies based on these operations.

At the same time, the well-foundedness of LTY implies well-foundedness results for BLP under the relevant operations. Thus the order type of BLP provides a lower bound for the strength of LTY. Once the BLP analysis reaches $\varepsilon_0$, the unprovability in PA follows naturally.

For me, this was an especially satisfying interaction between googology and mainstream mathematical logic: methods motivated by the analysis of extremely large ordinal notations produced a new way to attack a long-standing problem about Laver tables.

BLP also defines an ordinal notation system in its own right. Its structure is quite unlike the notation systems I had previously studied; I do not think I would have arrived at it without the Laver-table problem itself.

## What the theorem says

The main theorem of my paper is that the Laver-table function $F$ eventually dominates every computable function whose totality is provable in PA. In particular, the totality of $F$ is not provable in PA.

The paper also gives much stronger lower bounds for specific values associated with the fourth critical point. Some of these lower bounds already exceed Graham's number by a very large margin.

I believe the broader picture may be considerably stronger than the theorem currently proves. LTY appears to be an exceptionally strong explicit recursive ordinal notation system, and quantities associated with the Laver-table function---for example values often denoted by expressions such as $q(5)$ in the googology literature---appear to be vastly larger than Graham's number or $\text{TREE}(3)$. Claims at this scale should be understood carefully: some are rigorous lower bounds, while others depend on ongoing ordinal analysis and should presently be regarded as conjectural or heuristic.

Similarly, BLP may be among the strongest explicit ordinal notation systems for which substantial initial segments can still be expanded and analyzed by hand. The apparent well-foundedness of such concrete systems is also, at most, heuristic evidence concerning the consistency of the large-cardinal principles from which they arise; it is of course not a consistency proof.

The point is not that the current construction is the final word. BLP is still far from an optimal lower-bound analysis of the Laver-table function, and LTY appears to contain much more structure that we do not yet understand.

The AI-for-Math team at Peking University---Zhiyuan Zhang, Jiedong Jiang, and Leheng Chen---formalized the main result of my paper in Lean 4, providing an independent machine-checked verification of the central theorem.

---

# The rest of the story: ordinal analysis in the googology community

The remainder of this article concerns exploratory work in googology rather than the main theorem of the paper.

In googology, to **analyze** an ordinal notation usually means to build, often by a semi-enumerative process, a table matching expressions in a new notation with expressions in a better-understood system that are believed to denote the same ordinals. This procedure is not by itself a rigorous proof. However, in practice it is often the stage from which one discovers the precise embedding algorithm and the statements that can later be proved formally. Research on the $\Sigma_2$ levels of PoR, for example, visibly contains this pattern of experimentation followed by abstraction.

Because of an earlier experience involving questions of priority and attribution, I was initially cautious about discussing BLP publicly. I had previously released an outline of my proof of the well-foundedness of BMS and later answered detailed questions about the proof on Discord. A person who had challenged parts of my argument subsequently announced their own proof of BMS well-foundedness, and, in my view, its central idea was essentially the same as the one I had explained. After that experience, I preferred to discuss the early BLP work only with a small group of people I trusted, including [HypCos](https://github.com/hypcos), [24414-X357](https://www.zhihu.com/people/36bcfecf838dd64ab77957f5f54f4f4c), and [ProjectCF](https://www.zhihu.com/people/30430c133c10124a7d6c37ceb23e76d8).

As those discussions continued, I improved my BLP analyzer. At one point ProjectCF suggested that the `modify` operation might become easier to analyze if a single expansion were replaced by an infinite expansion, although no concrete rule was proposed at the time.

Later, under the writing guidance of my supervisor Noam Greenberg, I wrote the paper and [posted it on arXiv](https://arxiv.org/abs/2501.06733). After several revisions I felt comfortable discussing BLP publicly.

HypCos often analyzes notation systems using his web tool [Notation Explorer](https://hypcos.github.io/notation-explorer/). The tool requires a comparison algorithm for standard expressions---roughly, the expressions obtained by repeatedly expanding the maximal expression. I had not initially found such an algorithm, so automated analysis was blocked. Using my own expander, however, I analyzed BLP up to

$$
\operatorname{PTO}\bigl(\mathrm{KP}+\text{there exist }\Pi_n\text{-reflecting ordinals for all }n<\omega\bigr),
$$

corresponding to the BMS expression `(0,0,0)(1,1,1)(2,2,0)`.

Eventually I found the missing comparison algorithm. HypCos then began his own analysis and encountered the same difficulty with the `modify` operation. After asking me many detailed questions about the design principles behind BLP, he constructed a concrete variant in which `modify` was replaced by infinite expansion.

I called this system **IBLP (infinite basic Laver patterns)**. Unlike BLP, IBLP no longer directly inherits the well-foundedness proof from LTY. On the other hand, some of its local behavior is substantially easier to analyze.

A striking feature is that, throughout the analysis, IBLP behaves very much as though it were well founded, and at many natural milestones the BLP and IBLP expressions are extremely similar. For example, at the proof-theoretic ordinal often written $\operatorname{PTO}(\mathrm{KPM})$, corresponding to the BMS expression

```text
(0,0,0)(1,1,1)(2,1,1)(3,1,1)(3,1,0)(4,2,0)
```

the IBLP and BLP representatives are closely related.

<!-- TODO: restore the original IBLP and BLP figures for PTO(KPM). -->

I have long suspected that the rows of IBLP correspond to some mathematical structure in the way that rows of BLP correspond to elementary embeddings, but no such structure is currently known.

HypCos analyzed IBLP up to the BMS limit and somewhat beyond it. It is widely conjectured in the googology community that the BMS limit corresponds to something on the scale of $\operatorname{PTO}(Z_2)$. If that interpretation is correct and the analogy with the Laver-table function continues, it suggests the possibility of independence phenomena far beyond PA. This remains speculative: the relevant identification and independence statements have not been proved.

## Bad patterns and revised rules

The definition of BLP contains configurations that cannot be copied in the desired way. I informally call them **bad patterns**. My original treatment of these cases was deliberately crude, because I did not know the optimal rule; this weakens the resulting notation system.

A googologist known as ddfg is particularly good at finding descending sequences or structural failures in proposed ordinal notations. He studied IBLP and found a number of bad patterns, including unexpectedly small ones. Because BLP and IBLP are locally very similar, an IBLP bad pattern often points directly to a corresponding issue that must be handled in BLP.

There are many variants of these patterns, and they quickly become complicated. Worse, on the proof side, preserving the well-foundedness argument for BLP requires genuinely new ideas. After an intensive period of experimentation, I developed a replacement for the old `r-list` mechanism that I call **mark completion**. This led to improved rules for both BLP and IBLP.

The enhanced BLP analyzer is available [here](https://github.com/QiRenrui/basic-laver-pattern/blob/main/basic-laver-pattern-enhanced.py). A definition of the revised IBLP can be found on this [googology wiki page](https://wiki.googology.top/index.php/iBLP), and it also appears as IBLP in [Notation Explorer Rewritten](https://smilelee-lyx.github.io/ne-rewritten/). The revised BLP still has a well-foundedness proof, although this strengthened version has not yet been incorporated into the paper.

[Sigmoid](https://github.com/hzyhhzy) wrote a program to search automatically for bad patterns. With computations by ddfg and Sigmoid, many bad patterns in the revised IBLP were found; they occur much higher than the first bad patterns in the original version. Eventually the community isolated a particular IBLP pattern `(1,0)1(2,1,0)1(3,2,1,0)2(4,3,2)1(5,4,3,2)2(6,5,4)1` below which no bad pattern is currently known. There is no proof that none exists. For this experimental version of IBLP, I therefore defined the limit by truncating the system at that pattern.

One especially difficult family, which I call **stacked-row bad patterns**, remains unresolved. The smallest currently known example was again found by ddfg. There are still many basic structural questions open here.

## Beyond BMS: comparison with the Y-sequence

I continued the analysis from the point where HypCos stopped, using the **Y-sequence**, a well-known notation system in the googology community for ordinals beyond the usual BMS range. Early in this process, ddfg and [@投影序数](https://www.zhihu.com/people/8e40e4804405605bd10bdfbefc4898c6) provided substantial help.

My own analysis tended to make large jumps, so slower and more detailed independent calculations were important checks. Several milestones in my analysis were later independently confirmed, including

```text
Y(1,3,4,2,5,7,5)
Y(1,3,4,2,5,7,10,5)
Y(1,3,4,2,5,8)
Y(1,3,4,2,5,8,10)
Y(1,3,4,3)
Y(1,3,5)
Y(1,3,6)
Y(1,3,7)
Y(1,3,7,15)
Y(1,3,8)
Y(1,3,9)
Y(1,4)
```

as well as the limit of the Y-sequence.

[SmileLee (Li Yixiao)](https://github.com/SmileLee-lyx), an IMO 2018 gold medalist and a full-score gold medalist in the 2023 Alibaba Global Mathematics Competition, produced the strongest detailed comparison between IBLP and the Y-sequence. His analysis consists of 3,082 lines of explicit correspondences and essentially confirmed my earlier prediction that the Y-sequence limit corresponds in IBLP to

```text
(1,0)1(2,1,0)1(3,2,1,0)2(4,3,2)1(5,4,3,2)2
```

[@投影序数](https://www.zhihu.com/people/8e40e4804405605bd10bdfbefc4898c6) also contributed extensive discussion and later produced an extremely detailed independent analysis, including a fine analysis of the notoriously complicated node `Y(1,3,4,3)`.

[油手就行](https://www.zhihu.com/people/169dc1c7c92e544d3ab07438e3fbd7e5) independently carried out a detailed analysis through `Y(1,3,4,7)`.

At present, IBLP-related systems appear to be among the strongest ordinal notations for which the googology community has managed a substantial explicit analysis. They are also unusual in that they are not merely modifications of BMS or the Y-sequence, but arise from a different source: the combinatorics of Laver tables and elementary embeddings.

Even so, a fast-growing value such as $f_{\mathrm{IBLP}}(1000)$ may still be far smaller than the Laver-table quantities represented by values such as $q(5)$. There is a pleasing irony here: after years of developing ever stronger ordinal notation systems, googology may find that some of the enormous numbers it encountered much earlier through Laver tables lie far beyond the later constructions.

---

## Links

- Renrui Qi, [*Notes on Laver Tables*](https://arxiv.org/abs/2501.06733)
- Randall Dougherty, [*Critical points in an algebra of elementary embeddings*](https://arxiv.org/abs/math/9205202)
- [Enhanced basic Laver pattern analyzer](https://github.com/QiRenrui/basic-laver-pattern/blob/main/basic-laver-pattern-enhanced.py)
- [Notation Explorer rewritten](https://smilelee-lyx.github.io/ne-rewritten/)
- [IBLP wiki page](https://wiki.googology.top/index.php/iBLP)
