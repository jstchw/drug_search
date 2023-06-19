import { Container, Navbar } from 'react-bootstrap';
import { NavLink } from 'react-router-dom';
import HeaderLink from './HeaderLink'
import { Search as SearchIcon } from 'react-bootstrap-icons';

const Header = (props) => {
    return (
        // Navbar containing the links to all the pages of the webapp
        <Navbar bg="dark" variant="dark" sticky={'top'} className="header-full-width">
            <Container>
                <Navbar.Brand as={NavLink} to="/" onClick={e => window.location.href = '/'}>
                    <SearchIcon />DrugSearch
                </Navbar.Brand>
                <HeaderLink to="/dashboard">Dashboard</HeaderLink>
                <HeaderLink to="/settings">Settings</HeaderLink>
            </Container>
        </Navbar>
    );
}

export default Header;