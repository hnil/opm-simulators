-- The chapters link to each other as files, "02-examples.md#anchor", which
-- works when the markdown is browsed. In the single document the build makes,
-- the anchor alone is the target.
function Link(el)
  local anchor = el.target:match("^%d%d%-[%w%-]+%.md(#.+)$")
  if anchor then
    el.target = anchor
  end
  return el
end
