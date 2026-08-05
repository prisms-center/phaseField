#include <prismspf/utilities/assert.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <string>
#include <vector>

namespace Tags
{
  enum class Discretization
  {
    Explicit,
    Implicit
  };

  enum class Physics
  {
  };

  enum class PhaseField
  {
    MultiPhase,
  };

  enum class PDE
  {
    AllenCahn,
    CahnHilliard,
    Diffusion,
    Poisson,
    LinearElasticity
  };

  enum class Numerics
  {
    AMR,
    Nucleation,
    LinearSolve,
    NonlinearSolve,
    TimeDependentBC,
    Integration,
    GrainRemapping
  };
} // namespace Tags

/**
 * Application documentation struct
 */
struct ApplicationNode
{
  std::string           name;
  std::string           pretty_name;
  std::filesystem::path path;
  std::filesystem::path markdown;
  bool                  has_subfolders = false;
};

class ApplicationParser
{
public:
  explicit ApplicationParser(const std::string &file)
  {
    parse(file);
  }

  [[nodiscard]] const std::vector<ApplicationNode> &
  get_applications() const
  {
    return _applications;
  }

private:
  void
  parse(const std::filesystem::path &root)
  {
    if (!std::filesystem::exists(root) || !std::filesystem::is_directory(root))
      {
        return;
      }

    for (const auto &entry : std::filesystem::recursive_directory_iterator(root))
      {
        if (!entry.is_directory())
          {
            continue;
          }

        const auto app_dir = entry.path();

        if (!std::filesystem::is_regular_file(app_dir / "main.cc"))
          {
            continue;
          }

        ApplicationNode app;
        app.path = app_dir;

        for (const auto &sub : std::filesystem::directory_iterator(app_dir))
          {
            if (sub.is_regular_file() && sub.path().extension() == ".md")
              {
                app.markdown = sub.path();
                break;
              }
          }

        auto relative = std::filesystem::relative(app_dir, root);

        for (const auto &part : relative)
          {
            if (!app.name.empty())
              {
                app.name += "_";
              }

            app.name += part.string();
          }

        // Detect whether this application contains nested applications
        for (const auto &sub : std::filesystem::directory_iterator(app_dir))
          {
            if (sub.is_directory() && std::filesystem::exists(sub.path() / "main.cc"))
              {
                app.has_subfolders = true;
                break;
              }
          }

        _applications.push_back(std::move(app));
      }
  }

  std::vector<ApplicationNode> _applications;
};

class HTMLGenerator
{
public:
  explicit HTMLGenerator(const std::vector<ApplicationNode> &applications)
    : _applications(applications)
  {}

  void
  write(const std::string &input_filepath, const std::string &output_filepath) const
  {
    // Read the input file
    std::ifstream in(input_filepath);
    ASSERT(in.is_open(), "Could not open file", input_filepath);

    std::ostringstream buffer;
    buffer << in.rdbuf();
    std::string content = buffer.str();
    in.close();

    // Find the marker
    const std::string list_marker     = "@application_list";
    const std::size_t list_marker_pos = content.find(list_marker);
    ASSERT(list_marker_pos != std::string::npos,
           "Could not find list_marker in file",
           list_marker,
           input_filepath);

    // Replace the list_marker with the generated table
    content.replace(list_marker_pos, list_marker.size(), build_application_table());

    // Do it again for the application subpages
    const std::string subpage_marker     = "@application_subpages";
    const std::size_t subpage_marker_pos = content.find(subpage_marker);
    ASSERT(subpage_marker_pos != std::string::npos,
           "Could not find subpage_marker in file",
           subpage_marker,
           input_filepath);
    content.replace(subpage_marker_pos,
                    subpage_marker.size(),
                    build_application_subpages());

    // Write to output
    std::ofstream out(output_filepath);
    ASSERT(out.is_open(), "Could not open file", output_filepath);
    out << content;
    out.close();

    // Update the DoxygenLayout.xml
    // NOTE: This is hardcoded but I can't be bothered to change it
    const std::filesystem::path layout_input =
      std::filesystem::path(input_filepath).parent_path() / ".." / "DoxygenLayout.xml.in";

    const std::filesystem::path layout_output =
      std::filesystem::path(input_filepath).parent_path() / ".." / "DoxygenLayout.xml";

    update_layout(layout_input.string(), layout_output.string());
  }

private:
  void
  update_layout(const std::string &input_filepath,
                const std::string &output_filepath) const
  {
    std::ifstream in(input_filepath);
    ASSERT(in.is_open(), "Could not open file", input_filepath);

    std::ostringstream buffer;
    buffer << in.rdbuf();
    std::string content = buffer.str();
    in.close();

    const std::string marker = "@application_navigation";
    const auto        pos    = content.find(marker);
    ASSERT(pos != std::string::npos, "Could not find marker", marker, input_filepath);
    content.replace(pos, marker.size(), build_application_navigation());

    std::ofstream out(output_filepath);
    ASSERT(out.is_open(), "Could not open file", output_filepath);
    out << content;
    out.close();
  }

  [[nodiscard]] std::string
  build_application_table() const
  {
    std::ostringstream ss;

    ss << "<table>\n";
    ss << "  <thead>\n";
    ss << "    <tr>\n";
    ss << "      <th>Application</th>\n";
    ss << "    </tr>\n";
    ss << "  </thead>\n";
    ss << "  <tbody>\n";

    for (const auto &app : _applications)
      {
        ss << "    <tr>\n";
        ss << "      <td>\\ref " << app.name << " \"" << app.name << "\"</td>\n";
        ss << "    </tr>\n";
      }

    ss << "  </tbody>\n";
    ss << "</table>\n";

    return ss.str();
  }

  [[nodiscard]] std::string
  build_application_subpages() const
  {
    std::ostringstream ss;

    for (const auto &app : _applications)
      {
        ss << "/**\\page " << app.name << " " << app.name << "\n\n";

        std::ifstream in(app.markdown);
        ss << in.rdbuf();

        ss << "\n*/\n\n";
      }

    return ss.str();
  }

  [[nodiscard]] std::string
  build_application_navigation() const
  {
    std::ostringstream ss;

    for (const auto &app : _applications)
      {
        ss << "  <tab type=\"user\" "
           << "url=\"\\ref " << app.name << "\" "
           << "title=\"" << app.name << "\"/>\n";
      }

    return ss.str();
  }

  const std::vector<ApplicationNode> &_applications;
};

int
main(int argc, char *argv[])
{
  // Throw an error for the wrong script input
  if (argc != 4)
    {
      std::cerr
        << "Usage: \n  scrape_applications input_file output_file application_dir\n"
        << std::flush;
      std::exit(1);
    }

  // Check that the applications file exists
  const std::string in_file = argv[1];
  ASSERT(std::filesystem::exists(std::filesystem::path(in_file)),
         "Could not find file",
         in_file);
  const std::string out_file = argv[2];

  // Read applications into a more useful structure
  const std::string app_dir = argv[3];
  ApplicationParser parser(app_dir);
  const auto        applications = parser.get_applications();

  // Write to HTML
  HTMLGenerator writer(applications);
  writer.write(in_file, out_file);

  return 0;
}
