#set text(12pt, font: "New Computer Modern")
#set par(justify: true)
#show link: set text(fill: blue)

#text(style: "italic", size: 13pt)[Proceedings of the Royal Society A]

#align(center, text(18pt, weight: "bold")[
  Flight dynamics and optimal control of unpowered avian landings: Supplementary Material])

#align(center, text(13pt)[
  Tom Rottier\*, Phillipe Lavoie, Ben Parslew
])


#text(
  size: 10pt,
)[DOI: #link("https://doi.org/10.1098/rspa.2026-0103")[10.1098/rspa.2026-0103]]

#text(
  size: 10pt,
)[Corresponding author. Email: #link("mailto:tom.rottier@manchester.ac.uk")]


#let appendix(body) = {
  let app_state = state("")
  show heading.where(level: 1): it => {
    app_state.update(counter(heading).display())
    counter(math.equation).update(0)
    counter(figure.where(kind: image)).update(0)
    it
  }
  set heading(numbering: "A.1", supplement: [Supplementary Material])
  set math.equation(numbering: it => {
    [(#app_state.get()-#it)]
  })
  set figure(numbering: it => {
    [#app_state.get()-#it]
  })
  counter(heading).update(0)
  body
}
#show: appendix

#include "rigidbody.typ"
#pagebreak()

#include "opt-control.typ"
#pagebreak()

#bibliography(
  "../references.bib",
  style: "../proceedings-of-the-royal-society-b.csl",
)
