// Table layout and styling only; labels and data belong in slides.qmd.
#let presentation-table(columns, header-rows: 1, ..cells) = {
  set text(fill: rgb("22373b"), weight: "regular")
  show strong: set text(fill: rgb("22373b"), weight: "bold")
  show table.cell: cell => {
    set text(weight: "bold") if cell.y < header-rows
    cell
  }
  table(
    columns: columns,
    stroke: none,
    inset: (x: 6pt, y: 6pt),
    fill: none,
    table.hline(y: 0, stroke: 0.8pt + rgb("22373b")),
    table.hline(y: header-rows, stroke: 0.5pt + rgb("829497")),
    ..cells,
    table.hline(stroke: 0.8pt + rgb("22373b")),
  )
}

#let contingency-table(counts, row-label: none, column-label: none,
                       row-levels: (), column-levels: ()) = {
  assert(counts.len() == row-levels.len(), message: "One row label is required per data row.")
  assert(counts.all(row => row.len() == column-levels.len()),
         message: "One column label is required per data column.")
  let cells = ()
  for (i, row) in counts.enumerate() {
    if i == 0 {
      cells.push(table.cell(rowspan: counts.len(), align: center + horizon)[#row-label])
    }
    cells.push(table.cell(align: center)[#row-levels.at(i)])
    for value in row {
      cells.push(table.cell(align: center)[#value])
    }
  }
  presentation-table(
    (auto, auto, ..column-levels.map(_ => 1fr)),
    header-rows: 2,
    table.header(
      table.cell(colspan: 2, rowspan: 2)[],
      table.cell(colspan: column-levels.len(), align: center)[#column-label],
      ..column-levels.map(level => table.cell(align: center)[#level]),
    ),
    ..cells,
  )
}
