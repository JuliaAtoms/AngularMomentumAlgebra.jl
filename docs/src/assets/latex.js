requirejs.config({
    paths: {
        'mathjax': 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS_HTML',
    },
    shim: {
        'mathjax' : {
            exports: "MathJax"
        },
    }
});

require(['mathjax'], function(MathJax) {
    MathJax.Hub.Config({
        TeX: {
            Macros: {
                defd: "≝",
                abs: ["|#1|",1],
                ket: ["|#1\\rangle",1],
                bra: ["\\langle#1|",1],
                braket: ["\\langle#1|#2\\rangle",2],
                matrixel: ["\\langle#1|#2|#3\\rangle",3],
                redmatrixel: ["\\langle#1||#2||#3\\rangle",3],
                expect: ["\\langle#1\\rangle",1],
                vec: ["\\mathbf{#1}",1],
                mat: ["\\mathsf{#1}",1],
                conj: ["#1^*",1],
                im: "\\mathrm{i}",
                tensor : ["\\hat{\\mathbf{#1}}",1],
                operator : ["\\mathfrak{#1}",1],
                Hamiltonian : "\\operator{H}",
                hamiltonian : "\\operator{h}",
                Lagrangian : "\\operator{L}",
                fock : "\\operator{f}",
                lagrange : ["\\epsilon_{#1}",1],
                vary : ["\\delta_{#1}",1],
                onebody : ["(#1|#2)",2],
                twobody : ["[#1|#2]",2],
                twobodydx : ["[#1||#2]",2],
                direct : ["{\\operator{J}_{#1}}",1],
                exchange : ["{\\operator{K}_{#1}}",1],
                diff : ["\\mathrm{d}#1\\,",1],
                angroot : ["{\\prod}_{#1}", 1],
                wignerthreej : ["\\begin{pmatrix}#1\\end{pmatrix}", 1],
                wignersixj : ["\\begin{Bmatrix}#1\\end{Bmatrix}", 1],
                Heaviside: "\\Theta"
            }
        }
    });
})


