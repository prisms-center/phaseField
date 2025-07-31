#
# Try and find the deal.II library
#



set(DEAL_II_DIR "" CACHE PATH "An optional hint to a deal.II directory")
set_if_empty(DEAL_II_DIR "$ENV{DEAL_II_DIR}")

