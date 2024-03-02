import './Header.css';
import { Badge } from 'react-bootstrap';
import { Search as SearchIcon } from 'react-bootstrap-icons';
import { versionInfo } from '../../constants';

const Header = () => {
  return (
    <div className={'navbar-container p-0'}>
      <Badge className={'p-2'} style={{ cursor: 'pointer' }} onClick={() => (window.location.href = '/')}>
        <span className={'fs-5'}>
          <SearchIcon />
          {versionInfo.appName}
        </span>
        <Badge className={'version-sticker'}>{versionInfo.tag}</Badge>
      </Badge>
    </div>
  );
};

export default Header;
