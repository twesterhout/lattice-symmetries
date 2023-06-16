use LatticeSymmetries;

proc main() {
  var tree = new shared RoseTree(
    0,
    [ new shared RoseTree(1)
    , new shared RoseTree(2)
    , new shared RoseTree(
        3,
        [ new shared RoseTree(30)
        , new shared RoseTree(31)
        ]
      )
    ]
  );

  writeln(tree);
}
