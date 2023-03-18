import { NavLink } from 'react-router-dom';
import './HeaderLink.css'

const HeaderLink = ({ to, children }) => {
    return (
        <NavLink
            to={to}
            className="header-link"
        >
            {children}
        </NavLink>
    );
}

export default HeaderLink;