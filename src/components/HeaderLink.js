import { NavLink } from 'react-router-dom';
import Container from 'react-bootstrap/Container';

const HeaderLink = ({ to, children }) => {
    return (
        <Container>
            <NavLink to={to} className="header-link">{children}</NavLink>
        </Container>
    );
}

export default HeaderLink;