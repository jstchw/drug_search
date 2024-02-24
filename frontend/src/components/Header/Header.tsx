import "./Header.css";
import { Badge } from "react-bootstrap";
import { Search as SearchIcon } from "react-bootstrap-icons";

const Header = () => {
  return (
    <div className={"navbar-container p-0"}>
      <Badge
        className={"p-2"}
        style={{ cursor: "pointer" }}
        onClick={() => (window.location.href = "/")}
      >
        <span className={"fs-5"}>
          <SearchIcon />
          DrugSearch
        </span>
        <Badge className={"version-sticker"}>Alpha</Badge>
      </Badge>
    </div>
  );
};

export default Header;
