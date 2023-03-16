import { NavLink } from 'react-router-dom';
import './HeaderLink.css'

const HeaderLink = ({ to, children }) => {
    return (
        <NavLink
            to={to}
            exact
            className="header-link"
            activeClassName="header-link-active"
        >
            {children}
        </NavLink>
    );
}

export default HeaderLink;