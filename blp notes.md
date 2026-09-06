In iblp if you copy on $p,p+1,...,A_n[-L]$ instead of $p,p+1,...,n$, then you'll get the copying operation of blp with marks. In this way, it's easy to define blp with marks.



For a bls of a blp with marks, the only new requirement is that, if $i'\in M_i$, then let $tr_i(i')=(i_1,i_2,...,i_k)$, we have $j_i\equiv_{\theta_{i'+1}}j_{i_1}\circ j_{i_2}\circ ...\circ j_{i_{k-1}}$.

The main lemma of arxiv2501.96733 still holds. The only non-trivial thing we need to verify is definition 10, case (b). Let $tr_\sigma(y)=(i_1,i_2,...,i_l,i_{l+1},...,i_{l'})$, where $i_l=y_0$. We have $j_\sigma\equiv_{\theta_{y+1}}j_{i_1}\circ ... \circ j_{i_{l'-1}}$. If If $y_0<a_1$, then we have $j_n(j_{i_s})\equiv_{\theta_{a_1}}j_{i_s}$ for $s=l+1,l+2,...,l'-1$, applying $j_n$ to the expression we get the expression we need.

For another case, assume $tr_n(A_n[k+L])=(m_1,...,m_s)$, then $j_n\equiv_{\theta_{A_n[k+L]+1}}j_{m_1}\circ j_{m_2}\circ ...\circ j_{m_{s-1}}$. We also have $j_n(j_{i_{l+1}}\circ j_{i_{l+2}}\circ ... \circ j_{i_{l'-1}})\equiv_{j_n(j_{i_{l+1}}\circ j_{i_{l+2}}\circ ... \circ j_{i_{l'-1}})(\theta_{a_1})} j_n\circ j_{i_{l+1}}\circ j_{i_{l+2}}\circ ... \circ j_{i_{l'-1}}$.

Applying $j_n(j_{i_1}\circ ...\circ j_{i_l})$ on it, and note that $j_\tau(\theta_{a_1})\ge \theta_{x+1}$ by condition, we have $j_\tau\equiv_{\theta_{x+1}}j_n(j_{i_1})\circ j_n(j_{i_2})\circ ...\circ j_n(j_{i_l})\circ j_{m_1}\circ j_{m_2}\circ ...\circ j_{m_s-1}\circ j_{i_{l+1}}\circ j_{i_{l+2}}\circ ...\circ j_{i_{l'-1}}$,  which is what we need.



Ordinals of iblp is much easier to analyze than blp. We conjecture that it is well ordered, and some math object can interpret its rows instead of I3 embeddings.


