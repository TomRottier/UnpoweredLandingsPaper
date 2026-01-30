#set page(numbering: "1", margin: 0.5in)
#set heading(numbering: "1.")
#set text(12pt, font: "New Computer Modern")
#set math.equation(numbering: "(1)")
#set par(justify: true)

#align(center, text(18pt, weight: "bold")[
  Flight dynamics and optimal control of unpowered avian landings])

#align(center, text(13pt)[
  Tom Rottier, Phillipe Lavoie, Ben Parslew
])

#include "abstract.typ"

#include "introduction.typ"

#include "methods.typ"

= Results and Discussion
#include "results-dynamics.typ"
#include "results-opt-control.typ"
#include "results-constraints.typ"
#include "results-perching.typ"
#include "results-sensitivity.typ"

#include "conclusion.typ"


#pagebreak()

#bibliography("references.bib", style: "proceedings-of-the-royal-society-b.csl")

#pagebreak()

= Supplementary Material <supp-mat>
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
#include "appendix/rigidbody.typ"
#include "appendix/opt-control.typ"
// #include "appendix/energy-height.typ"

// == Supplementary material A <supp-mat-rigidbody>
// == Supplementary material B <supp-mat-optcontrol>
// == Supplementary material C <supp-mat-energyheight>

